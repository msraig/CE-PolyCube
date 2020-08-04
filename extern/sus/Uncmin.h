/*! \file Uncmin.h

  \brief
  C++ implementation of UNCMIN routines for unconstrained minimization
  
  This software is available at 
  http://www.smallwaters.com/software/cpp/uncmin.html
 
  Version 0080309
  This release updates original work by Brad Hanson (http://www.b-a-h.com/) in 2001.
 
  The UNCMIN routines are based on Appendix A of:
 
  Dennis, J. E., & Schnabel, R. B. (1996). Numerical methods for
  unconstrained optimization and nonlinear equations. Philadelphia:
  Society for Industrual and Applied Mathematics.
 
  The original UNCMIN routines in FORTRAN are available from GAMS
  (http://gams.nist.gov/) as OPTIF9 and OPTIF0.
 
  FORTRAN UNCMIN routines were translated to java (Uncmin_f77.java) by
  Steve Verrill (steve@ws13.fpl.fs.fed.us).
  Uncmin_F77.java is available at http://www1.fpl.fs.fed.us/optimization.html

  The routines in Uncmin_f77.java were translated to C++ by
  Brad Hanson.
  
  (c) 2008, Werner Wothke, maintenance programming and documentation.
 
 
  Uncmin is a templated class with template parameters:
 
  V - A vector class of double precision numbers. Member functions
      that should be defined for class V are:
 
      V& operator=(const V &A); // assignment
 
      V& operator=(const double& scalar); // assign a value to all elements
 
      double& operator()(int i); // FORTRAN style one-offset subscripting
 
      int size(); // returns number of elements in vector
 
  M - A matrix class of double precision numbers. Member functions that
      should be defined for class M are:
 
      double operator()(int i, int j); // FORTRAN style one-offset subscripting
 
      int num_rows(); // returns number of rows in matrix
      int num_cols(); // returns number of columns in matrix

  FUNC - A class defining the function to minimize, its gradient and hessian.
         An example of how this class should be defined is:
  
    // V is vector class, M is matrix class
    template <class V, class M>
    class minclass
    {
       double f_to_minimize(V &x);
       // Returns value of function to minimize
 
       void gradient(V &x, V &g);
       // Returns gradient of function in g. This does not have
       // to be defined, in which case the gradient is approximated.
       // If the gradient is not defined Analytic_Gradient should 
       // return 0;
 
       void hessian(V &x, M &h);
       // Returns hessian of function in h. This does not have
       // to be defined, in which case the hessian is approximated.
       // If the hessian is not defined Analytic_Hessian should 
       // return 0;
 
       int HasAnalyticGradient();
       // Returns 1 if gradient is defined, otherwise returns 0
 
       int HasAnalyticHessian();
       // Returns 1 if hessian is defined, otherwise returns 0
 
       int ValidParameters(V &x);
       // Returns 1 if f_to_minimize(x) is defined for argument x.
 
       int dim();
       // Returns dimension of problem (number of elements in vector x
       // passed to f_to_minimize)
    }
 
    An example of how to use Uncmin to minimize a function defined in
    a class minclass with vector and matrix classes from the Simple
    C++ Numerical Toolkit (http://www.smallwaters.com/software/cpp/scppnt.html):
 
    #include "Uncmin.h"
    #include "vec.h"  // Simple C++ Numeric Toolkit (SCPPNT) vector class
    #include "cmat.h" // Simple C++ Numeric Toolkit (SCPPNT) matrix class
    #include <cstdio>
 
    using namespace SCPPNT;
 
    // Dimension of problem (number of parameters)
    const int dim = 4;
 
    // Create function object
    minclass<Vector<double>, Matrix<double> > test;
 
    // create Uncmin object
    Uncmin<Vector<double>, Matrix<double>, minclass> min(&test);
 
    // Starting values
    double start_values[4] = {1.0, 1.0, 1.0, 1.0}; // use some appropriate values here
    Vector<double> start(start_values, start_values+dim);

    // xpls will contain solution, gpls will contain
    // gradient at solution after call to min.Minimize
    Vector<double> xpls(dim), gpls(dim);

    // fpls contains the value of the function at
    // the solution given by xpls.
    double fpls;
 
    // Minimize function.
    // Minimize returns zero if successful (0, 1, 2, 3
    // returned by GetMessage), non-zero if not successful.
    if (min.Minimize(start, xpls,  fpls,  gpls))
    {
       int msg = min.GetMessage(); // find out what went wrong
       std::printf("\nMessage returned from Uncmin: %d\n", msg);
    }

    The GetMessage member function returns the status of the solution found
    by the last call to the Minimize member function.
    Possible values returned by GetMessage are:
 
    0 = Optimal solution found, terminated with gradient small
    1 = Terminated with gradient small, xpls is probably optimal.
    2 = Terminated with stepsize small, xpls is probably optimal.
    3 = Lower point cannot be found, xpls is probably optimal.
    4 = Iteration limit (default 150) exceeded.
    5 = Too many large steps, function may be unbounded.
    -1 = Analytic gradient check requested, but no analytic gradient supplied
    -2 = Analytic hessian check requested, but no analytic hessian supplied
    -3 = Illegal dimension
    -4 = Illegal tolerance
    -5 = Illegal iteration limit
    -6 = Minimization function has no good digits
    -7 = Iteration limit exceeded in lnsrch
    -20 = Function not defined at starting value
    -21 = Check of analytic gradient failed
    -22 = Check of analytic hessian failed


    The symbol BOOST_NO_STDC_NAMESPACE should be defined if the C++ compiler
    does not put the standard C library functions in the std namespace
    (this is needed for Visual C++ 6).
 
    The symbol BOOST_NO_LIMITS should be defined if the standard C++
    library header file <limits> is not available.
    (this is needed for Visual C++ 6).

 
    Version History
 
    Version 080309 (March 9, 2008)
    
    Edited doxygen-compatible documentation, second pass. Added descriptive comments, 
    as much as possible, and incorporated any useful information from the existing 
    FORTRAN and Java documentation sets. Removed these out-dated documentations 
    afterwards.
    
    Version 071226 (December 26, 2007)
 
    Refactored num_cols property to num_columns in UNCMIN::minimize().
    Tagged all public methods and properties with doxygen-compatible comments.
 
    Version 0.001206 (February 10, 2001)
 
    Added BOOST_NO_STDC_NAMESPACE and BOOST_NO_LIMITS symbols which
    need to be defined to compile Uncmin.h using Microsoft Visual C++ 6.

 */

#ifndef UNCMIN_H_
#define UNCMIN_H_

#include <cstdio>
#include <cmath>

#define USING_PROJ 1

#ifdef BOOST_NO_STDC_NAMESPACE
// If the standard C library functions are not in the std namespace (like Microsoft Visual C++ 6)
namespace std
{ using ::sqrt; using ::pow; using ::fabs; using ::log; using ::FILE; using ::fprintf;}
#endif

#ifndef BOOST_NO_LIMITS
// include limits for epsilon function used in SetTolerances.
#include <limits>
#else
// If there is no limits header.
#include <float.h>
#endif


//CGAL part

#include <Eigen/Dense>
#include <vector>
#include <list>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_triangle_primitive.h>
#include <CGAL/Polyhedron_incremental_builder_3.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Polyhedron_items_with_id_3.h>
#include <CGAL/AABB_face_graph_triangle_primitive.h>

using Eigen::Vector3d;
typedef CGAL::Simple_cartesian<double> Ker;
typedef Ker::FT FT;
typedef Ker::Point_3 Point3;
typedef Ker::Point_2 Point2;
typedef std::vector<FT> Scalar_vector;
//typedef CGAL::Barycentric_coordinates::Triangle_coordinates_2<Ker> Triangle_coordinates;

typedef Ker::Segment_3 Segment;
typedef CGAL::Polyhedron_3<Ker, CGAL::Polyhedron_items_with_id_3> Polyhedron;
typedef CGAL::AABB_face_graph_triangle_primitive<Polyhedron> Primitive;
typedef CGAL::AABB_traits<Ker, Primitive> Traits;
typedef CGAL::AABB_tree<Traits> Tree;
typedef Tree::Point_and_primitive_id Point_and_primitive_id;
typedef Polyhedron::HalfedgeDS             HalfedgeDS;

Tree *search_tree = NULL;
Polyhedron polyhedron;
std::vector<Vector3d> proj_normals;

//std::vector<Vector3d> proj_normals;


/*! Documentation of member functions is given in their definitions which follow
 the Uncmin class definition */
template<class V, class M, class FUNC> class Uncmin
{

public:

  typedef FUNC function_type;   //!< Type definition of function object to minimize.

  /* Constructor. */
  Uncmin(FUNC *f);
  
  /* Constructor. */
  Uncmin(int d = 1);

  /* Destructor */
  ~Uncmin();

  /* Calculate function minimum. */
  int Minimize(V &start, V &xpls, double &fpls, V &gpls, M &hpls);

  int Minimize(V &start, V &xpls, double &fpls, V &gpls, M &hpls, int normal_x, int normal_y, int normal_z);

  /* Version of Minimize where hessian matrix is allocated locally rather
   than passed as an argument to the function. */
  int Minimize(V &start, V &xpls, double &fpls, V &gpls);

  /* Get object containing function to minimize, gradient, and hessian. */
  FUNC *GetFunction();

  /* Set object containing function to minimize, gradient, and hessian. */
  void SetFunction(FUNC *f);

  /* Methods to set options. */
  void SetFuncExpensive(int iexp);

  /* Set typical magnitude of each argument to function. */
  int SetScaleArg(V typsiz);

  /* Set typical magnitude of function near minimum. */
  void SetScaleFunc(double fscale);

  /* Set the number of reliable digits returned by f_to_minimize. */
  void SetNumDigits(double ndigit);

  /* Set maximum number of iterations. */
  void SetMaxIter(int maxiter);

  /* Set step tolerance, gradient tolerance, and machine epsilon. */
  void SetTolerances(double step, double grad, double macheps = -1.0);

  /* Set maximum allowable scaled step length at any iteration. */
  void SetMaxStep(double stepmx);

  /* Set initial trust region radius. */
  void SetTrustRegionRadius(double dlt);

  /* Set flags to check analytic gradient and hessian. */
  void SetChecks(int check_gradient, int check_hessian = 0);

  /* Set flags to indicate whether results are printed. */
  void SetPrint(std::FILE *file, int print_results, int print_iterations = 0);

  /* Set method to use to solve minimization problem.
   1 = line search, 2 = double dogleg, 3 = More-Hebdon.  */
  void SetMethod(int method);

  /* Return stopping criteria computed at the end of the last call to Minimize. */
  void GetStoppingCriteria(double &gradtol, double &steptol, int &iterations);

  /*! Returns message generated in the last call to Minimize. */
  int GetMessage();

private:

  FUNC *minclass; //!<  Class that contains function to miminize, gradient, and hessian. 


  int algorithm; //!<  Indicates algorithm to use to solve minimization problem: 1 = line search, 2 = double dogleg, 3 = More-Hebdon.

  std::FILE *mfile; //!< File to print messages to.

  int n; //!< Dimension of problem, provided in minclass by minclass->dim().

  double epsm; //!< Machine epsilon.


  double mLastGradCrit;  //!< Stopping criterion for gradient.
  

  double mLastStepCrit; //!<  Stopping criterion for argument change.

  
  int mLastIteration; //!<  Stopping criterion by number of iterations.


  int mLastMessage; //!<  Message generated by last call to Minimize.

  /* Options */

  int iexp;  //!< Flag indicating computational complexity. Set iexp = 1, for expensive calculations, iexp = 0 otherwise. If iexp = 1, the Hessian will be evaluated by secant (BFGS) update.
  
  /* These two flags are set from minclass using minclass.Analytic_Gradient() and
   minclass.Analytic_Methods() */
  int iagflg; //!< Flag: iagflag = 0 if an analytic gradient is not supplied, provided by minclass.
  int iahflg; //!< Flag: iahflag = 0 if an analytic Hessian is not supplied, provided by minclass.

  int fCheckGradient; //!< fCheckGradient = 0 if analytic gradient is not checked.

  int fCheckHessian; //!< fCheckHessian = 0 if analytic hessian is not checked.

  int fPrintResults; //!< fPrintResults = 0 if results are not printed to mfile.

  int fPrintIterationResults; //!< fPrintIterationResults = 0 if results of each iteration are not printed to mfile.

  /* Tolerances */

  V typsiz; //!< Typical size of of each argument of function to minimize.

  double fscale; //!< Estimate of scale of objective function.

  int ndigit; //!< Number of good digits in the minimization function. Special case, set ndigit = -1, for function providing within 1 or 2 of full number of significant digits.

  int itnlim; //!< Maximum number of allowable iterations.

  double mdlt; //!< Trust region radius. Special case,  mdlt = -1, uses length of initial scaled gradient instead.

  double gradtl; //!< Tolerance at which the gradient is considered close enough to zero to terminate the algorithm.

  double steptl; //!< Tolerance at which scaled distance between two successive iterates is considered close enough to zero to terminate algorithm. 
  
  double stepmx; //!< Maximum allowable step size.

  /* Private functions used in minimization. */

  void dfault();

  int optchk(V &x, V &sx, int &msg, int &method);

  void fstofd(V &xpls, double fpls, V &g, V &sx, double rnoise);

  void fstofd(V &xpls, double fpls, V &g, V &sx, double rnoise, int normal_x, int normal_y, int normal_z);

  void fstofd(V &xpls, V &fpls, M &a, V &sx, double rnoise, V &fhat);

  int grdchk(V &x, double f, V &g, V &typsiz, V &sx, double fscale, double rnf, double analtl,
      V &gest, int &msg);

  void optstp(V &xpls, double fpls, V &gpls, V &x, int &itncnt, int &icscmx, int &itrmcd, V &sx,
      double &fscale, int &itnlim, int &iretcd, bool &mxtake, int &msg);

  void hsnint(M &a, V &sx, int method);

  void sndofd(V &xpls, double &fpls, M &a, V &sx, double &rnoise, V &stepsz, V &anbr);

  void sndofd(V &xpls, double &fpls, M &a, V &sx, double &rnoise, V &stepsz, V &anbr, int normal_x, int normal_y, int normal_z);

  int heschk(V &x, double &f, V &g, M &a, V &typsiz, V &sx, double rnf, double analtl, int &iagflg,
      V &udiag, V &wrk1, V &wrk2, int &msg);

  int heschk(V &x, double &f, V &g, M &a, V &typsiz, V &sx, double rnf, double analtl, int &iagflg,
	  V &udiag, V &wrk1, V &wrk2, int &msg, int normal_x, int normal_y, int normal_z);

  void result(V &x, double &f, V &g, M &a, V &p, int &itncnt, int iflg);

  void bakslv(M &a, V &x, V &b);

  void chlhsn(M &a, V &sx, V &udiag);

  void choldc(M &a, double diagmx, double tol, double &addmax);

  void dogdrv(V &x, double &f, V &g, M &a, V &p, V &xpls, double &fpls, V &sx, double &dlt,
      int &iretcd, bool &mxtake, V &sc, V &wrk1, V &wrk2, V &wrk3);

  void dogstp(V &g, M &a, V &p, V &sx, double rnwtln, double &dlt, bool &nwtake, bool &fstdog,
      V &ssd, V &v, double &cln, double &eta, V &sc);

  void forslv(M &a, V &x, V &b);

  void fstocd(V &x, V &sx, double rnoise, V &g);

  void hookdr(V &x, double &f, V &g, M &a, V &udiag, V &p, V &xpls, double &fpls, V &sx,
      double &dlt, int &iretcd, bool &mxtake, double &amu, double &dltp, double &phi,
      double &phip0, V &sc, V &xplsp, V &wrk0, int &itncnt);

  void hookst(V &g, M &a, V &udiag, V &p, V &sx, double rnwtln, double &dlt, double &amu,
      double &dltp, double &phi, double &phip0, bool &fstime, V &sc, bool &nwtake, V &wrk0);

  void lltslv(M &a, V &x, V &b);

  void lnsrch(V &x, double &f, V &g, V &p, V &xpls, double &fpls, bool &mxtake, int &iretcd, V &sx);

  void lnsrch(V &x, double &f, V &g, V &p, V &xpls, double &fpls, bool &mxtake, int &iretcd, V &sx, int normal_x, int normal_y, int normal_z);


  void mvmltl(M &a, V &x, V &y);

  void mvmlts(M &a, V &x, V &y);

  void mvmltu(M &a, V &x, V &y);

  void qraux1(M &r, int i);

  void qraux2(M &r, int i, double a, double b);

  void qrupdt(M &a, V &u, V &v);

  void secfac(V &x, V &g, M &a, V &xpls, V &gpls, int &itncnt, double rnf, int &iagflg,
      bool &noupdt, V &s, V &y, V &u, V &w);

  void secunf(V &x, V &g, M &a, V &udiag, V &xpls, V &gpls, int &itncnt, double rnf, int &iagflg,
      bool &noupdt, V &s, V &y, V &t);

  void tregup(V &x, double &f, V &g, M &a, V &sc, V &sx, bool &nwtake, double &dlt, int &iretcd,
      V &xplsp, double &fplsp, V &xpls, double &fpls, bool &mxtake, int method, V &udiag);

  double ddot(V &x, V &y);

  double dnrm2(V &x);

  double max(double a, double b);

  double min(double a, double b);

};

/*!
  \brief 
  Constructor that sets dimension and function to minimize.
  
  The minclass function object must define:
                    
  1. a method, f_to_minimize, to minimize. f_to_minimize must have the form 
         public static double f_to_minimize(double x[])
     where x is the vector of arguments to the function and the return value is the 
     value of the function evaluated at x.  

  2. a method, gradient, that has the form
         public static void gradient(double x[], double g[])
     where g is the gradient of f evaluated at x.  This method will have an empty body if 
     the user does not wish to provide an analytic estimate of the gradient.
     
  3. a method, hessian, that has the form
         public static void hessian(double x[],double h[][])
     where h is the Hessian of f evaluated at x.  This method will have an empty body if the 
     user does not wish to provide an analytic estimate of the Hessian. If the user wants 
     Uncmin++ to check the Hessian, then the hessian method should only fill the lower 
     triangle (and diagonal)of h.
  
   Minimization parameters are set at default values.
  
  \section template_args Template Parameters
  
  \param M  SCPPNT Matrix type.
  \param V  SCPPNT Vector type.
  \param FUNC Type of function to minimize.
  
  \section function_args Function Parameters
   
  \param[in]  *f  Pointer to function object f.
 */
template<class V, class M, class FUNC> Uncmin<V, M, FUNC>::Uncmin(FUNC *f) :
  minclass(f)
{
  if (minclass)
    n = minclass->dim();
  else
    n = 0;

  mfile = 0;

  mLastGradCrit = 0.0;
  mLastStepCrit = 0.0;
  mLastIteration = 0;
  mLastMessage = 0;

  /* Set default parameter values */
  dfault();
}

/*!
  \brief 
  Constructor that sets dimension, but not function to minimize.
  
  Note: Function to minimize must be set by calling SetFunction. 
  Minimization parameters are set to default values.
  
  \section template_args Template Parameters
  
  \param M  SCPPNT Matrix type.
  \param V  SCPPNT Vector type.
  \param FUNC Type of function to minimize.
  
  \section function_args Function Parameters
   
  \param[in]  d  Dimension d (number of function arguments).
 */
template<class V, class M, class FUNC> Uncmin<V, M, FUNC>::Uncmin(int d) :
  minclass(0), typsiz(d)
{

  n = d;

  mfile = 0;

  mLastGradCrit = 0.0;
  mLastStepCrit = 0.0;
  mLastIteration = 0;
  mLastMessage = 0;

  /* Set default parameter values */
  dfault();

}

/*!
	\brief Destructor of uncmin.
 */
template<class V, class M, class FUNC> Uncmin<V, M, FUNC>::~Uncmin()
{

}

/*!
  \brief 
  Assigns flag indicating whether function is expensive to evaluate.
  
  Set iexp =1 if function is expensive to evaluate, = 0 otherwise.  
  If iexp = 1, then the Hessian will be evaluated by secant (BFGS) update 
  rather than analytically or by finite differences.

  Default value is 0 if this function is not called.
  
  \section template_args Template Parameters
  
  \param M  SCPPNT Matrix type.
  \param V  SCPPNT Vector type.
  \param FUNC Type of function to minimize.
  
  \section function_args Function Parameters
   
  \param[in]  iexp  Set iexp =1 if function is expensive to evaluate, iexp = 0 otherwise.
 */
template<class V, class M, class FUNC> void Uncmin<V, M, FUNC>::SetFuncExpensive(int iexp)
{
  this->iexp = iexp;
}

/*!
  \brief
  Assigns typical magnitude of each argument to function.

  Returns zero if dimension matches, otherwise returns 1.
 
  The default is to set all elements of typsiz to 1 if this function is not called.
  
  \section template_args Template Parameters
  
  \param M  SCPPNT Matrix type.
  \param V  SCPPNT Vector type.
  \param FUNC Type of function to minimize.
  
  \section function_args Function Parameters
   
  \param[in]  typsiz  Vector describing the magnitude of all function arguments.
 */
template<class V, class M, class FUNC> int Uncmin<V, M, FUNC>::SetScaleArg(V typsiz)
{
  if (typsiz.size() != n)
    return 1;

  this->typsiz = typsiz;

  return 0;
}

/*!
  \brief
  Assigns typical magnitude of function near minimum.
  
  \section template_args Template Parameters
  
  \param M  SCPPNT Matrix type.
  \param V  SCPPNT Vector type.
  \param FUNC Type of function to minimize.
  
  \section function_args Function Parameters
   
  \param[in]  fscale  Magnitude of function values near minimum.
                      A default value of 1 is used for fscale if this
                      function is not called.
 */
template<class V, class M, class FUNC> void Uncmin<V, M, FUNC>::SetScaleFunc(double fscale)
{
  this->fscale = fscale;
}

/*!
  \brief
  Assigns the number of reliable digits returned by f_to_minimize.
  
  If the argument to the function is -1 this
  means that f_to_minimize is expected to
  provide within one or two of the full number of
  significant digits available.
 
  The default value of -1 is used for ndigit when
  this function is not called.

  \section template_args Template Parameters
  
  \param M  SCPPNT Matrix type.
  \param V  SCPPNT Vector type.
  \param FUNC Type of function to minimize.
  
  \section function_args Function Parameters
   
  \param[in]  ndigit  Number of digits precision required for the solution.
 */
template<class V, class M, class FUNC> void Uncmin<V, M, FUNC>::SetNumDigits(double ndigit)
{
  this->ndigit = ndigit;
}

/*!
  \brief
  Assigns the maximum number of iterations.
  
  The default value of 150 is used if this function is not called.

  \section template_args Template Parameters
  
  \param M  SCPPNT Matrix type.
  \param V  SCPPNT Vector type.
  \param FUNC Type of function to minimize.
  
  \section function_args Function Parameters
   
  \param[in]  maxiter  Maximum number of iterations.
 */
template<class V, class M, class FUNC> void Uncmin<V, M, FUNC>::SetMaxIter(int maxiter)
{
  this->itnlim = maxiter;
}

/*!
  \brief
  Assigns step tolerance, gradient tolerance, and the machine epsilon.
  
  Note: For any of these values <= zero, default values are used.
  
  \section template_args Template Parameters
  
  \param M  SCPPNT Matrix type.
  \param V  SCPPNT Vector type.
  \param FUNC Type of function to minimize.
  
  \section function_args Function Parameters
   
 	\param[in]  step  Tolerance at which scaled distance between two successive iterates
   	is considered close enough to zero to terminate algorithm. If the value of
  	step is <= zero a default value is used.
 	\param[in]  grad  Tolerance at which gradient is considered close enough
  	to zero to terminate algorithm. If the value of
  	grad is <= zero a default value is used.
  \param[in]  macheps Largest relative spacing. Machine epsilon == 2**(-T+1) where
  	(1 + 2**-T) == 1. This argument has a default value of -1 if not included
  	in function call.
 */
template<class V, class M, class FUNC> void Uncmin<V, M, FUNC>::SetTolerances(double step,
    double grad, double macheps)
{
  if (macheps <= 0.0) // use default value
  {

#ifdef BOOST_NO_LIMITS
    epsm = DBL_EPSILON; // from float.h
#else
    // use epsilon from numeric_limits in standard C++ library
    epsm = std::numeric_limits<double>::epsilon();
#endif

  }
  else
  {
    epsm = macheps;
  }

  if (grad <= 0.0)
  {
    gradtl = std::pow(epsm, 1.0/3.0); // default value
  }
  else
  {
    gradtl = grad;
  }

  if (step <= 0.0)
  {
    steptl = std::sqrt(epsm); // default value
  }
  else
  {
    steptl = step;
  }

}

/*!
  \brief
  Assigns maximum allowable scaled step length at any iteration.
  
  The maximum step length is used to prevent steps that would cause the 
  optimization algorithm to overflow or leave the domain of iterest,
  as well as to detect divergence. It should be chosen small enough
  to prevent the first two of these occurrences, but larger
  than any anticipated reasonable step size. The algorithm will
  halt if it takes steps larger than the maximum allowable step
  length on 5 consecutive iterations.
 
  If the argument to the function is <= 0 then a default value of
  stepmx will be computed in optchk.
 
  A default value of 0 is used for stepmx when this function is
  not called.

  \section template_args Template Parameters
  
  \param M  SCPPNT Matrix type.
  \param V  SCPPNT Vector type.
  \param FUNC Type of function to minimize.
  
  \section function_args Function Parameters
   
  \param[in]  stepmx  Maximum (scaled) step size.
 */
template<class V, class M, class FUNC> void Uncmin<V, M, FUNC>::SetMaxStep(double stepmx)
{
  this->stepmx = stepmx;
}

/*!
  \brief
  Assigns initial trust region radius.
  
  Used when method = 2 or 3, ignored when method = 1. 
  The value should be what the user considers a reasonable 
  scaled step length for the first iteration, and should be 
  less than the maximum step length.
 
  If dlt = -1, then the length of the initial scaled gradient is used instead.
 
  A default value of -1 is used if this function is not called.
 
  \section template_args Template Parameters
  
  \param M  SCPPNT Matrix type.
  \param V  SCPPNT Vector type.
  \param FUNC Type of function to minimize.
  
  \section function_args Function Parameters
   
  \param[in]  dlt  Maximum (scaled) step size of first iteration.
 */
template<class V, class M, class FUNC> void Uncmin<V, M, FUNC>::SetTrustRegionRadius(double dlt)
{
  this->mdlt = dlt;
}

/*!
  \brief
  Returns pointer to object containing function to minimize, gradient, and hessian.
 
  \section template_args Template Parameters
  
  \param M  SCPPNT Matrix type.
  \param V  SCPPNT Vector type.
  \param FUNC Type of function to minimize.
 */
template<class V, class M, class FUNC> FUNC * Uncmin<V, M, FUNC>::GetFunction()
{
  return minclass;
}


/*!
  \brief
  Assigns pointer to object containing function to minimize, gradient, and hessian.

  \section template_args Template Parameters
  
  \param M  SCPPNT Matrix type.
  \param V  SCPPNT Vector type.
  \param FUNC Type of function to minimize.
  
  \section function_args Function Parameters
   
  \param[in]  *f  Pointer to function object f.
 */
template<class V, class M, class FUNC> void Uncmin<V, M, FUNC>::SetFunction(FUNC *f)
{
  minclass = f;

  /* If dimension does not match change n and re-intialize typsiz */
  if (n != f->dim())
  {
    n = f->dim();
    V defaultTypSiz(n, 1.0);
    typsiz = defaultTypSiz;
  }
}

/*!
  \brief
  Assigns flags to check analytic gradient and hessian.

  The default is to not check the gradient and hessian, if this function is not called
  (fCheckGradient = 0 and fCheckHessian = 0).
 
  Note that check_hessian can be left off function call in which case default value of
  0 is used.
  
  \section template_args Template Parameters
  
  \param M  SCPPNT Matrix type.
  \param V  SCPPNT Vector type.
  \param FUNC Type of function to minimize.
  
  \section function_args Function Parameters
   
  \param[in]  check_gradient  Flag to check analytic gradient of function: 1 = yes, 0 = no.
  \param[in]  check_hessian  Flag to check analytic hessian of function: 1 = yes, 0 = no.
 */
template<class V, class M, class FUNC> void Uncmin<V, M, FUNC>::SetChecks(int check_gradient,
    int check_hessian)
{
  fCheckGradient = check_gradient;
  fCheckHessian = check_hessian;
}

/*!
  \brief
  Assigns flags to indicate whether results are printed and specify which file is 
  receive the print output.

  Printing is done to open file indicated by first argument.

  The default action is not to print any output.
  
  \section template_args Template Parameters
  
  \param M  SCPPNT Matrix type.
  \param V  SCPPNT Vector type.
  \param FUNC Type of function to minimize.
  \param FILE Type of file for writing print output to.
  
  \section function_args Function Parameters
   
  \param[in]  *file  Pointer to the open file that will receive the printed results.
  \param[in]  print_results  Flag (if non-zero) to print standard results.
  \param[in]  print_iterations  Flag (if non-zero) to print intermediate results at each 
      iteration (this argument can be left off function call in which case the
      default value of 0 is used).
  */
template<class V, class M, class FUNC> void Uncmin<V, M, FUNC>::SetPrint(std::FILE *file,
    int print_results, int print_iterations)
{

  mfile = file;

  if (!file)
  {
    fPrintResults = 0;
    fPrintIterationResults = 0;
  }
  else
  {
    fPrintResults = print_results;
    fPrintIterationResults = print_iterations;
  }

}

/*!
  \brief
  Assigns method to use to solve minimization problem.
  
  A default value is 1 (line search) if SetMethod is not called.

  \section template_args Template Parameters
  
  \param M  SCPPNT Matrix type.
  \param V  SCPPNT Vector type.
  \param FUNC Type of function to minimize.
  
  \section function_args Function Parameters
   
  \param[in]  method  Method indicator: 1 = line search, 2 = double dogleg, 3 = More-Hebdon.
 */
template<class V, class M, class FUNC> void Uncmin<V, M, FUNC>::SetMethod(int method)
{
  if (method == 2 || method == 3)
    algorithm = method;
  else
    algorithm = 1; // use default value of 1 for any values of argument other than 2 or 3
}

/*!
  \brief
  Returns stopping criteria computed at the end of the last call to Minimize.

  \section template_args Template Parameters
  
  \param M  SCPPNT Matrix type.
  \param V  SCPPNT Vector type.
  \param FUNC Type of function to minimize.
  
  \section function_args Function Parameters
   
  \param[out]  &gradtol  Address of gradient tolerance, i.e., how close the gradient is to zero.
  \param[out]  &steptol  Address of maximum scaled difference between function parameters in 
      consecutive iterations.
  \param[out]  &iterations  Address of number of iterations used.
 */
template<class V, class M, class FUNC> void Uncmin<V, M, FUNC>::GetStoppingCriteria(
    double &gradtol, double &steptol, int &iterations)
{
  gradtol = mLastGradCrit;
  steptol = mLastStepCrit;
  iterations = mLastIteration;
}

/*!
  \brief
  Returns message index generated in the last call to Minimize.
 
  Possible values of message are:
 
    0 = Optimal solution found, terminated with gradient small.
    1 = Terminated with gradient small, xpls is probably optimal.
    2 = Terminated with stepsize small, xpls is probably optimal.
    3 = Lower point cannot be found, xpls is probably optimal.
    4 = Iteration limit (default 150) exceeded.
    5 = Too many large steps, function may be unbounded.
   -1 = Analytic gradient check requested, but no analytic gradient supplied
   -2 = Analytic hessian check requested, but no analytic hessian supplied
   -3 = Illegal dimension
   -4 = Illegal tolerance
   -5 = Illegal iteration limit
   -6 = Minimization function has no good digits
   -7 = Iteration limit exceeded in lnsrch
  -20 = Function not defined at starting value
  -21 = Check of analytic gradient failed
  -22 = Check of analytic hessian failed

  \section template_args Template Parameters
  
  \param M  SCPPNT Matrix type.
  \param V  SCPPNT Vector type.
  \param FUNC Type of function to minimize.
 */
template<class V, class M, class FUNC> inline int Uncmin<V, M, FUNC>::GetMessage()
{
  return mLastMessage;
}

/*!
  \brief
  Method to compute values that minimize the function.

  Translated by Brad Hanson from optdrv_f77 by Steve Verrill.
  
  Returns 0 if mLastMessage is 0, 1, 2, or 3 (indicating optimal solution probably found), 
  otherwise returns non-zero.
  
  On exit data member mLastMessage is set to indicate status of solution.

  Possible values of mLastMessage are:
    0 = Optimal solution found, terminated with gradient small.
    1 = Terminated with gradient small, xpls is probably optimal.
    2 = Terminated with stepsize small, xpls is probably optimal.
    3 = Lower point cannot be found, xpls is probably optimal.
    4 = Iteration limit (default 150) exceeded.
    5 = Too many large steps, function may be unbounded.
   -1 = Analytic gradient check requested, but no analytic gradient supplied.
   -2 = Analytic hessian check requested, but no analytic hessian supplied.
   -3 = Illegal dimension.
   -4 = Illegal tolerance.
   -5 = Illegal iteration limit.
   -6 = Minimization function has no good digits. 
   -7 = Iteration limit exceeded in lnsrch.
  -20 = Function not defined at starting value.
  -21 = Check of analytic gradient failed.
  -22 = Check of analytic hessian failed.

  \section template_args Template Parameters
  
  \param M  SCPPNT Matrix type.
  \param V  SCPPNT Vector type.
  \param FUNC Type of function to minimize.
  
  \section function_args Function Parameters
   
  \param[in]   &start  Address of vector containing initial estimate of location of minimum.
  \param[out]  &xpls  Address of vector containing local minimum (size must be n on input).
  \param[out]  &fpls  Address of function value at local minimum.
  \param[out]  &gpls  Address of vector containing gradient at local minimum (size must be 
      n on input). 
  \param[out]  &hpls  Address of matrix containing hessian at local minimum (size must be 
      n X n on input). Only the lower half should be input; the upper off-diagonal entries 
      can be set to zero. On return, the lower half of hpls contains the Cholesky factor 
      of the Hessian, evaluated at the solution (i.e., H = LL').
 */
template<class V, class M, class FUNC> int Uncmin<V, M, FUNC>::Minimize(V &start, V &xpls,
    double &fpls, V &gpls, M &hpls)
{

  bool noupdt;
  bool mxtake;

  int i, j, num5, remain, ilow, ihigh;
  int icscmx;
  int iretcd;
  int itncnt;

  /*
   	Return code.
   		itrmcd =  0:    Optimal solution found
   		itrmcd =  1:    Terminated with gradient small, X is probably optimal
   		itrmcd =  2:    Terminated with stepsize small, X is probably optimal
   		itrmcd =  3:    Lower point cannot be found, X is probably optimal
   		itrmcd =  4:    Iteration limit (150) exceeded
   		itrmcd =  5:    Too many large steps, function may be unbounded.
   */
  int itrmcd;

  double rnf, analtl, dltsav, amusav, dlpsav, phisav, phpsav;

  double f;

  double amu = 0.0;
  double dltp = 0.0;
  double phi = 0.0;
  double phip0 = 0.0;

  int method = algorithm;

  mLastGradCrit = 0.0;
  mLastStepCrit = 0.0;
  mLastIteration = 0;
  mLastMessage = -9999;

  iagflg = minclass->HasAnalyticGradient();
  iahflg = minclass->HasAnalyticHessian();

  /*
    msg
   	Flag to set various options.
   	Bits which indicate flags which control checks of gradient and hessian and what output is printed.
   	msg should be zero or a sum of the following values which indicate which flags to set or not set.

   		1 = Turn off check of analytic gradient.
   		2 = Turn off check of analytic hessian.
   		4 = Turn off printing of results.
   		8 = Turn on printing of results at each iteration.
   */
  int msg = (fCheckGradient == 0) * 1 + (fCheckHessian == 0) * 2 + (fPrintResults == 0) * 4
      + (fPrintIterationResults != 0) * 8;

  /* Check that dimensions match */
  if (start.size() != n || xpls.size() != n || gpls.size() != n || hpls.num_rows() != n
      || hpls.num_columns() != n) // Refactored function "num_columns", ww, 12-22-2007
  {
    mLastMessage = -3;
    return -1;
  }

  /* Check that starting values are valid */
  if (!(minclass->ValidParameters(start)))
  {
    mLastMessage = -20;
    return -1;
  }

  /* Workspace */
  M &a = hpls;
  V udiag(n); //workspace (for diagonal of Hessian)
  V g(n); // workspace (for gradient at current iterate)
  V p(n); // workspace for step
  V sx(n); // workspace (for scaling vector)
  V wrk0(n), wrk1(n), wrk2(n), wrk3(n);

  V x(start);

  double dlt = mdlt;

  dltsav = amusav = dlpsav = phisav = phpsav = 0.0;

  // INITIALIZATION

  for (i = 1; i <= n; i++)
  {
    p(i) = 0.0;
  }

  itncnt = 0;
  iretcd = -1;

  /* Check for valid options */
  if (optchk(x, sx, msg, method))
  {
    mLastMessage = msg;
    return -1;
  }

  rnf = max(std::pow(10.0, -ndigit), epsm);
  analtl = max(.01, std::sqrt(rnf));

  if (!(msg & 4))
  {
    num5 = n/5;
    remain = n%5;

    std::fprintf(mfile, "\n\nOPTDRV          Typical x\n\n");

    ilow = -4;
    ihigh = 0;

    for (i = 1; i <= num5; i++)
    {
      ilow += 5;
      ihigh += 5;

      std::fprintf(mfile, "%d--%d     ", ilow, ihigh);

      for (j = 1; j <= 5; j++)
      {
        std::fprintf(mfile, "%lf  ", (double) typsiz(ilow+j-1));
      }
      std::fprintf(mfile, "\n");
    }

    ilow += 5;
    ihigh = ilow + remain - 1;

    std::fprintf(mfile, "%d--%d     ", ilow, ihigh);

    for (j = 1; j <= remain; j++)
    {
      std::fprintf(mfile, "%lf  ", (double) typsiz(ilow+j-1));
    }

    std::fprintf(mfile, "\n");

    std::fprintf(mfile, "\n\nOPTDRV      Scaling vector for x\n\n");

    ilow = -4;
    ihigh = 0;

    for (i = 1; i <= num5; i++)
    {
      ilow += 5;
      ihigh += 5;

      std::fprintf(mfile, "%d--%d     ", ilow, ihigh);

      for (j = 1; j <= 5; j++)
      {
        std::fprintf(mfile, "%lf  ", (double) sx(ilow+j-1));
      }
      std::fprintf(mfile, "\n");
    }

    ilow += 5;
    ihigh = ilow + remain - 1;

    std::fprintf(mfile, "%d--%d     ", ilow, ihigh);

    for (j = 1; j <= remain; j++)
    {
      std::fprintf(mfile, "%lf  ", (double) sx(ilow+j-1));
    }

    std::fprintf(mfile, "\n");

    std::fprintf(mfile, "\n\nOPTDRV      Typical f = %f\n", fscale);
    std::fprintf(mfile, "OPTDRV      Number of good digits in");
    std::fprintf(mfile, " f_to_minimize = %d\n", ndigit);
    std::fprintf(mfile, "OPTDRV      Gradient flag");
    std::fprintf(mfile, " = %d\n", iagflg);
    std::fprintf(mfile, "OPTDRV      Hessian flag");
    std::fprintf(mfile, " = %d\n", iahflg);
    std::fprintf(mfile, "OPTDRV      Expensive function calculation flag");
    std::fprintf(mfile, " = %d\n", iexp);
    std::fprintf(mfile, "OPTDRV      Method to use");
    std::fprintf(mfile, " = %d\n", method);
    std::fprintf(mfile, "OPTDRV      Iteration limit");
    std::fprintf(mfile, " = %d\n", itnlim);
    std::fprintf(mfile, "OPTDRV      Machine epsilon");
    std::fprintf(mfile, " = %e\n", epsm);
    std::fprintf(mfile, "OPTDRV      Maximum step size");
    std::fprintf(mfile, " = %f\n", stepmx);
    std::fprintf(mfile, "OPTDRV      Step tolerance");
    std::fprintf(mfile, " = %e\n", steptl);
    std::fprintf(mfile, "OPTDRV      Gradient tolerance");
    std::fprintf(mfile, " = %e\n", gradtl);
    std::fprintf(mfile, "OPTDRV      Trust region radius");
    std::fprintf(mfile, " = %f\n", dlt);
    std::fprintf(mfile, "OPTDRV      Relative noise in");
    std::fprintf(mfile, " f_to_minimize = %e\n", rnf);
    std::fprintf(mfile, "OPTDRV      Analytical fd tolerance");
    std::fprintf(mfile, " = %e\n", analtl);
  }

  // EVALUATE FCN(X)
  f = minclass->f_to_minimize(x);

  // EVALUATE ANALYTIC OR FINITE DIFFERENCE GRADIENT AND CHECK ANALYTIC
  // GRADIENT, IF REQUESTED.
  if (iagflg == 0)
  {
    fstofd(x, f, g, sx, rnf);
	 // fstofd(x, f, g, sx, rnf, n);
  }
  else
  {
    minclass->gradient(x, g);

    if (!(msg & 1))
    {
      if (grdchk(x, f, g, typsiz, sx, fscale, rnf, analtl, wrk1, msg))
      {
        mLastMessage = msg;
        return -1; // msg = -21
      }
    }
  }

  optstp(x, f, g, wrk1, itncnt, icscmx, itrmcd, sx, fscale, itnlim, iretcd, mxtake, msg);

  if (itrmcd != 0)
  {
    fpls = f;

    for (i = 1; i <= n; i++)
    {
      xpls(i) = x(i);
      gpls(i) = g(i);
    }
  }
  else
  {
    if (iexp == 1)
    {
      // IF OPTIMIZATION FUNCTION EXPENSIVE TO EVALUATE (IEXP=1), THEN
      // HESSIAN WILL BE OBTAINED BY SECANT (BFGS) UPDATES.  GET INITIAL HESSIAN.

      hsnint(a, sx, method);
    }
    else
    {
      // EVALUATE ANALYTIC OR FINITE DIFFERENCE HESSIAN AND CHECK ANALYTIC
      // HESSIAN IF REQUESTED (ONLY IF USER-SUPPLIED ANALYTIC HESSIAN
      // ROUTINE minclass->hessian FILLS ONLY LOWER TRIANGULAR PART AND DIAGONAL OF A).

      if (iahflg == 0)
      {
        if (iagflg == 1)
        {
          fstofd(x, g, a, sx, rnf, wrk1);
        }
        else
        {
          sndofd(x, f, a, sx, rnf, wrk1, wrk2);
        }
      }
      else
      {
        if (msg & 2)
        {
          minclass->hessian(x, a);
        }
        else
        {
          if (heschk(x, f, g, a, typsiz, sx, rnf, analtl, iagflg, udiag, wrk1, wrk2, msg))
          {
            mLastMessage = msg;
            return -1; // msg = -22
          }
          // HESCHK EVALUATES minclass->hessian AND CHECKS IT AGAINST THE FINITE
          // DIFFERENCE HESSIAN WHICH IT CALCULATES BY CALLING FSTOFD
          // (IF IAGFLG .EQ. 1) OR SNDOFD (OTHERWISE).
        }
      }
    }

    if (!(msg & 4))
      result(x, f, g, a, p, itncnt, 1);

    // ITERATION
    while (itrmcd == 0)
    {
      itncnt++;
      mLastIteration = itncnt; // BAH

      // FIND PERTURBED LOCAL MODEL HESSIAN AND ITS LL+ DECOMPOSITION
      // (SKIP THIS STEP IF LINE SEARCH OR DOGSTEP TECHNIQUES BEING USED WITH
      // SECANT UPDATES.  CHOLESKY DECOMPOSITION L ALREADY OBTAINED FROM
      // SECFAC.)
      if (iexp != 1 || method == 3)
      {
        chlhsn(a, sx, udiag);
      }

      // SOLVE FOR NEWTON STEP:  AP=-G
      for (i = 1; i <= n; i++)
      {
        wrk1(i) = -g(i);
      }

      lltslv(a, p, wrk1);

      // DECIDE WHETHER TO ACCEPT NEWTON STEP  XPLS=X + P
      // OR TO CHOOSE XPLS BY A GLOBAL STRATEGY.
      if (iagflg == 0 && method != 1)
      {
        dltsav = dlt;

        if (method != 2)
        {
          amusav = amu;
          dlpsav = dltp;
          phisav = phi;
          phpsav = phip0;
        }
      }

      if (method == 1)
      {
        lnsrch(x, f, g, p, xpls, fpls, mxtake, iretcd, sx);
        if (iretcd == -7) // BAH
        {
          mLastMessage = iretcd;
          return -1;
        }
      }
      else if (method == 2)
      {
        dogdrv(x, f, g, a, p, xpls, fpls, sx, dlt, iretcd, mxtake, wrk0, wrk1, wrk2, wrk3);
      }
      else
      {
        hookdr(x, f, g, a, udiag, p, xpls, fpls, sx, dlt, iretcd, mxtake, amu, dltp, phi, phip0,
            wrk0, wrk1, wrk2, itncnt);
      }

      // IF COULD NOT FIND SATISFACTORY STEP AND FORWARD DIFFERENCE
      // GRADIENT WAS USED, RETRY USING CENTRAL DIFFERENCE GRADIENT.
      if (iretcd == 1 && iagflg == 0)
      {
        // SET IAGFLG FOR CENTRAL DIFFERENCES
        iagflg = -1;

        if (!(msg & 4))
        {
          std::fprintf(mfile, "\nOPTDRV      Shift from forward to central");
          std::fprintf(mfile, " differences");
          std::fprintf(mfile, "\nOPTDRV      in iteration %d", itncnt);
          std::fprintf(mfile, "\n");
        }

        fstocd(x, sx, rnf, g);

        if (method == 1)
        {
          // SOLVE FOR NEWTON STEP:  AP=-G
          for (i = 1; i <= n; i++)
          {
            wrk1(i) = -g(i);
          }

          lltslv(a, p, wrk1);

          lnsrch(x, f, g, p, xpls, fpls, mxtake, iretcd, sx);
          if (iretcd == -7) // BAH
          {
            mLastMessage = iretcd;
            return -1;
          }
        }
        else
        {
          dlt = dltsav;
          if (method == 2)
          {
            // SOLVE FOR NEWTON STEP:  AP=-G
            for (i = 1; i <= n; i++)
            {
              wrk1(i) = -g(i);
            }

            lltslv(a, p, wrk1);

            dogdrv(x, f, g, a, p, xpls, fpls, sx, dlt, iretcd, mxtake, wrk0, wrk1, wrk2, wrk3);
          }
          else
          {
            amu = amusav;
            dltp = dlpsav;
            phi = phisav;
            phip0 = phpsav;

            chlhsn(a, sx, udiag);

            // SOLVE FOR NEWTON STEP:  AP=-G
            for (i = 1; i <= n; i++)
            {
              wrk1(i) = -g(i);
            }

            lltslv(a, p, wrk1);

            hookdr(x, f, g, a, udiag, p, xpls, fpls, sx, dlt, iretcd, mxtake, amu, dltp, phi,
                phip0, wrk0, wrk1, wrk2, itncnt);
          }
        }
      }
      // CALCULATE STEP FOR OUTPUT
      for (i = 1; i <= n; i++)
      {
        p(i) = xpls(i) - x(i);
      }

      // CALCULATE GRADIENT AT XPLS
      if (iagflg == -1)
      {
        // CENTRAL DIFFERENCE GRADIENT
        fstocd(xpls, sx, rnf, gpls);
      }
      else if (iagflg == 0)
      {
        // FORWARD DIFFERENCE GRADIENT
        fstofd(xpls, fpls, gpls, sx, rnf);
      }
      else
      {
        // ANALYTIC GRADIENT
        minclass->gradient(xpls, gpls);
      }

      // CHECK WHETHER STOPPING CRITERIA SATISFIED
      optstp(xpls, fpls, gpls, x, itncnt, icscmx, itrmcd, sx, fscale, itnlim, iretcd, mxtake, msg);

      if (itrmcd == 0)
      {
        // EVALUATE HESSIAN AT XPLS
        if (iexp != 0)
        {
          if (method == 3)
          {
            secunf(x, g, a, udiag, xpls, gpls, itncnt, rnf, iagflg, noupdt, wrk1, wrk2, wrk3);
          }
          else
          {
            secfac(x, g, a, xpls, gpls, itncnt, rnf, iagflg, noupdt, wrk0, wrk1, wrk2, wrk3);
          }
        }
        else
        {
          if (iahflg == 1)
          {
            minclass->hessian(xpls, a);
          }
          else
          {
            if (iagflg == 1)
            {
              fstofd(xpls, gpls, a, sx, rnf, wrk1);
            }
            else
            {
              sndofd(xpls, fpls, a, sx, rnf, wrk1, wrk2);
            }
          }
        }

        if (msg & 8)
        {
          result(xpls, fpls, gpls, a, p, itncnt, 1);
        }

        // X <-- XPLS  AND  G <-- GPLS  AND  F <-- FPLS
        f = fpls;

        for (i = 1; i <= n; i++)
        {
          x(i) = xpls(i);
          g(i) = gpls(i);
        }
      }
    }

    // TERMINATION
    // -----------
    // RESET XPLS,FPLS,GPLS,  IF PREVIOUS ITERATE SOLUTION
    //
    if (itrmcd == 3)
    {
      fpls = f;

      for (i = 1; i <= n; i++)
      {
        xpls(i) = x(i);
        gpls(i) = g(i);
      }
    }
  }

  // PRINT RESULTS
  if (!(msg & 4))
  {
    result(xpls, fpls, gpls, a, p, itncnt, 0);
  }
  mLastMessage = itrmcd;
  return (itrmcd >= 0 && itrmcd <= 3) ? 0 : 1;
}


template<class V, class M, class FUNC> int Uncmin<V, M, FUNC>::Minimize(V &start, V &xpls,
	double &fpls, V &gpls, M &hpls, int normal_x, int normal_y, int normal_z)
{

	bool noupdt;
	bool mxtake;

	int i, j, num5, remain, ilow, ihigh;
	int icscmx;
	int iretcd;
	int itncnt;

	/*
	  Return code.
		  itrmcd =  0:    Optimal solution found
		  itrmcd =  1:    Terminated with gradient small, X is probably optimal
		  itrmcd =  2:    Terminated with stepsize small, X is probably optimal
		  itrmcd =  3:    Lower point cannot be found, X is probably optimal
		  itrmcd =  4:    Iteration limit (150) exceeded
		  itrmcd =  5:    Too many large steps, function may be unbounded.
	 */
	int itrmcd;

	double rnf, analtl, dltsav, amusav, dlpsav, phisav, phpsav;

	double f;

	double amu = 0.0;
	double dltp = 0.0;
	double phi = 0.0;
	double phip0 = 0.0;

	int method = algorithm;

	mLastGradCrit = 0.0;
	mLastStepCrit = 0.0;
	mLastIteration = 0;
	mLastMessage = -9999;

	//both 1
	iagflg = minclass->HasAnalyticGradient();
	iahflg = minclass->HasAnalyticHessian();

	/*
	  msg
	  Flag to set various options.
	  Bits which indicate flags which control checks of gradient and hessian and what output is printed.
	  msg should be zero or a sum of the following values which indicate which flags to set or not set.

		  1 = Turn off check of analytic gradient.
		  2 = Turn off check of analytic hessian.
		  4 = Turn off printing of results.
		  8 = Turn on printing of results at each iteration.
	 */
	int msg = (fCheckGradient == 0) * 1 + (fCheckHessian == 0) * 2 + (fPrintResults == 0) * 4
		+ (fPrintIterationResults != 0) * 8;

	/* Check that dimensions match */
	if (start.size() != n || xpls.size() != n || gpls.size() != n || hpls.num_rows() != n
		|| hpls.num_columns() != n) // Refactored function "num_columns", ww, 12-22-2007
	{
		mLastMessage = -3;
		return -1;
	}

	/* Check that starting values are valid */
	if (!(minclass->ValidParameters(start)))
	{
		mLastMessage = -20;
		return -1;
	}

	/* Workspace */
	M &a = hpls;
	V udiag(n); //workspace (for diagonal of Hessian)
	V g(n); // workspace (for gradient at current iterate)
	V p(n); // workspace for step
	V sx(n); // workspace (for scaling vector)
	V wrk0(n), wrk1(n), wrk2(n), wrk3(n);

	V x(start);

	double dlt = mdlt;

	dltsav = amusav = dlpsav = phisav = phpsav = 0.0;

	// INITIALIZATION

	for (i = 1; i <= n; i++)
	{
		p(i) = 0.0;
	}

	itncnt = 0;
	iretcd = -1;

	/* Check for valid options */
	if (optchk(x, sx, msg, method))
	{
		mLastMessage = msg;
		return -1;
	}

	rnf = max(std::pow(10.0, -ndigit), epsm);
	analtl = max(.01, std::sqrt(rnf));

	if (!(msg & 4))
	{
		num5 = n / 5;
		remain = n % 5;

		std::fprintf(mfile, "\n\nOPTDRV          Typical x\n\n");

		ilow = -4;
		ihigh = 0;

		for (i = 1; i <= num5; i++)
		{
			ilow += 5;
			ihigh += 5;

			std::fprintf(mfile, "%d--%d     ", ilow, ihigh);

			for (j = 1; j <= 5; j++)
			{
				std::fprintf(mfile, "%lf  ", (double)typsiz(ilow + j - 1));
			}
			std::fprintf(mfile, "\n");
		}

		ilow += 5;
		ihigh = ilow + remain - 1;

		std::fprintf(mfile, "%d--%d     ", ilow, ihigh);

		for (j = 1; j <= remain; j++)
		{
			std::fprintf(mfile, "%lf  ", (double)typsiz(ilow + j - 1));
		}

		std::fprintf(mfile, "\n");

		std::fprintf(mfile, "\n\nOPTDRV      Scaling vector for x\n\n");

		ilow = -4;
		ihigh = 0;

		for (i = 1; i <= num5; i++)
		{
			ilow += 5;
			ihigh += 5;

			std::fprintf(mfile, "%d--%d     ", ilow, ihigh);

			for (j = 1; j <= 5; j++)
			{
				std::fprintf(mfile, "%lf  ", (double)sx(ilow + j - 1));
			}
			std::fprintf(mfile, "\n");
		}

		ilow += 5;
		ihigh = ilow + remain - 1;

		std::fprintf(mfile, "%d--%d     ", ilow, ihigh);

		for (j = 1; j <= remain; j++)
		{
			std::fprintf(mfile, "%lf  ", (double)sx(ilow + j - 1));
		}

		std::fprintf(mfile, "\n");

		std::fprintf(mfile, "\n\nOPTDRV      Typical f = %f\n", fscale);
		std::fprintf(mfile, "OPTDRV      Number of good digits in");
		std::fprintf(mfile, " f_to_minimize = %d\n", ndigit);
		std::fprintf(mfile, "OPTDRV      Gradient flag");
		std::fprintf(mfile, " = %d\n", iagflg);
		std::fprintf(mfile, "OPTDRV      Hessian flag");
		std::fprintf(mfile, " = %d\n", iahflg);
		std::fprintf(mfile, "OPTDRV      Expensive function calculation flag");
		std::fprintf(mfile, " = %d\n", iexp);
		std::fprintf(mfile, "OPTDRV      Method to use");
		std::fprintf(mfile, " = %d\n", method);
		std::fprintf(mfile, "OPTDRV      Iteration limit");
		std::fprintf(mfile, " = %d\n", itnlim);
		std::fprintf(mfile, "OPTDRV      Machine epsilon");
		std::fprintf(mfile, " = %e\n", epsm);
		std::fprintf(mfile, "OPTDRV      Maximum step size");
		std::fprintf(mfile, " = %f\n", stepmx);
		std::fprintf(mfile, "OPTDRV      Step tolerance");
		std::fprintf(mfile, " = %e\n", steptl);
		std::fprintf(mfile, "OPTDRV      Gradient tolerance");
		std::fprintf(mfile, " = %e\n", gradtl);
		std::fprintf(mfile, "OPTDRV      Trust region radius");
		std::fprintf(mfile, " = %f\n", dlt);
		std::fprintf(mfile, "OPTDRV      Relative noise in");
		std::fprintf(mfile, " f_to_minimize = %e\n", rnf);
		std::fprintf(mfile, "OPTDRV      Analytical fd tolerance");
		std::fprintf(mfile, " = %e\n", analtl);
	}

	// EVALUATE FCN(X)
	//f = minclass->f_to_minimize(x);
	f = minclass->f_to_minimize(x, normal_x, normal_y, normal_z);

	//f = minclass->f_to_minimize(x, normal_x, normal_y, normal_z);

	// EVALUATE ANALYTIC OR FINITE DIFFERENCE GRADIENT AND CHECK ANALYTIC
	// GRADIENT, IF REQUESTED.
	if (iagflg == 0)
	{
		//fstofd(x, f, g, sx, rnf);
		fstofd(x, f, g, sx, rnf, normal_x, normal_y, normal_z);
	}
	else
	{
		minclass->gradient(x, g);

		if (!(msg & 1))
		{
			if (grdchk(x, f, g, typsiz, sx, fscale, rnf, analtl, wrk1, msg))
			{
				mLastMessage = msg;
				return -1; // msg = -21
			}
		}
	}

	optstp(x, f, g, wrk1, itncnt, icscmx, itrmcd, sx, fscale, itnlim, iretcd, mxtake, msg);

	if (itrmcd != 0)
	{
		fpls = f;

		for (i = 1; i <= n; i++)
		{
			xpls(i) = x(i);
			gpls(i) = g(i);
		}
	}
	else
	{
		//itrmcd = 0:    Optimal solution found
		if (iexp == 1)
		{
			// IF OPTIMIZATION FUNCTION EXPENSIVE TO EVALUATE (IEXP=1), THEN
			// HESSIAN WILL BE OBTAINED BY SECANT (BFGS) UPDATES.  GET INITIAL HESSIAN.

			hsnint(a, sx, method);
		}
		else
		{
			// EVALUATE ANALYTIC OR FINITE DIFFERENCE HESSIAN AND CHECK ANALYTIC
			// HESSIAN IF REQUESTED (ONLY IF USER-SUPPLIED ANALYTIC HESSIAN
			// ROUTINE minclass->hessian FILLS ONLY LOWER TRIANGULAR PART AND DIAGONAL OF A).

			if (iahflg == 0)
			{
				if (iagflg == 1)
				{
					fstofd(x, g, a, sx, rnf, wrk1);
					//fstofd(x, g, a, sx, rnf, wrk1, normal_x, normal_y, normal_z);
				}
				else
				{
					//sndofd(x, f, a, sx, rnf, wrk1, wrk2);
					sndofd(x, f, a, sx, rnf, wrk1, wrk2, normal_x, normal_y, normal_z);
				}
			}
			else
			{
				if (msg & 2)
				{
					minclass->hessian(x, a);
				}
				else
				{
					//if (heschk(x, f, g, a, typsiz, sx, rnf, analtl, iagflg, udiag, wrk1, wrk2, msg))
					if (heschk(x, f, g, a, typsiz, sx, rnf, analtl, iagflg, udiag, wrk1, wrk2, msg, normal_x, normal_y, normal_z))
					{
						mLastMessage = msg;
						return -1; // msg = -22
					}
					// HESCHK EVALUATES minclass->hessian AND CHECKS IT AGAINST THE FINITE
					// DIFFERENCE HESSIAN WHICH IT CALCULATES BY CALLING FSTOFD
					// (IF IAGFLG .EQ. 1) OR SNDOFD (OTHERWISE).
				}
			}
		}

		if (!(msg & 4))
			result(x, f, g, a, p, itncnt, 1);

		// ITERATION
		while (itrmcd == 0)
		{
			itncnt++;
			mLastIteration = itncnt; // BAH

			// FIND PERTURBED LOCAL MODEL HESSIAN AND ITS LL+ DECOMPOSITION
			// (SKIP THIS STEP IF LINE SEARCH OR DOGSTEP TECHNIQUES BEING USED WITH
			// SECANT UPDATES.  CHOLESKY DECOMPOSITION L ALREADY OBTAINED FROM
			// SECFAC.)
			if (iexp != 1 || method == 3)
			{
				chlhsn(a, sx, udiag);
			}

			// SOLVE FOR NEWTON STEP:  AP=-G
			for (i = 1; i <= n; i++)
			{
				wrk1(i) = -g(i);
			}

			lltslv(a, p, wrk1);

			// DECIDE WHETHER TO ACCEPT NEWTON STEP  XPLS=X + P
			// OR TO CHOOSE XPLS BY A GLOBAL STRATEGY.
			if (iagflg == 0 && method != 1)
			{
				dltsav = dlt;

				if (method != 2)
				{
					amusav = amu;
					dlpsav = dltp;
					phisav = phi;
					phpsav = phip0;
				}
			}

			if (method == 1)
			{
				/*if (normal_x > 0)
				{
					p(0) = 0;
				}
				if (normal_y > 0)
				{
					p(1) = 0;
				}
				if (normal_z > 0)
				{
					p(2) = 0;
				}*/
#if USING_PROJ
				//main newton step
				/*if (normal_x > 0)
				{
					Point3 query(x(1), x(2), x(3));
					Point_and_primitive_id pp = search_tree->closest_point_and_primitive(query);
					int face_id = pp.second->id();

					Vector3d tmp_g(p(1), p(2), p(3));
					double l = proj_normals[face_id].dot(tmp_g);
					Vector3d final_g = tmp_g - l * proj_normals[face_id];
					p(1) = final_g(0);
					p(2) = final_g(1);
					p(3) = final_g(2);
				}*/


#endif
				



				//lnsrch(x, f, g, p, xpls, fpls, mxtake, iretcd, sx);
				lnsrch(x, f, g, p, xpls, fpls, mxtake, iretcd, sx, normal_x, normal_y, normal_z);
				if (iretcd == -7) // BAH
				{
					mLastMessage = iretcd;
					return -1;
				}
			}
			else if (method == 2)
			{
				dogdrv(x, f, g, a, p, xpls, fpls, sx, dlt, iretcd, mxtake, wrk0, wrk1, wrk2, wrk3);
			}
			else
			{
				hookdr(x, f, g, a, udiag, p, xpls, fpls, sx, dlt, iretcd, mxtake, amu, dltp, phi, phip0,
					wrk0, wrk1, wrk2, itncnt);
			}

			// IF COULD NOT FIND SATISFACTORY STEP AND FORWARD DIFFERENCE
			// GRADIENT WAS USED, RETRY USING CENTRAL DIFFERENCE GRADIENT.
			if (iretcd == 1 && iagflg == 0)
			{
				// SET IAGFLG FOR CENTRAL DIFFERENCES
				iagflg = -1;

				if (!(msg & 4))
				{
					std::fprintf(mfile, "\nOPTDRV      Shift from forward to central");
					std::fprintf(mfile, " differences");
					std::fprintf(mfile, "\nOPTDRV      in iteration %d", itncnt);
					std::fprintf(mfile, "\n");
				}

				fstocd(x, sx, rnf, g);

				if (method == 1)
				{
					// SOLVE FOR NEWTON STEP:  AP=-G
					for (i = 1; i <= n; i++)
					{
						wrk1(i) = -g(i);
					}

					lltslv(a, p, wrk1);

					/*if (normal_x > 0)
					{
						p(0) = 0;
					}
					if (normal_y > 0)
					{
						p(1) = 0;
					}
					if (normal_z > 0)
					{
						p(2) = 0;
					}*/

					//try changing newton direction here
					//starting from 1
#if USING_PROJ
					//
					/*if (normal_x > 0)
					{
						Point3 query(x(1), x(2), x(3));
						Point_and_primitive_id pp = search_tree->closest_point_and_primitive(query);
						int face_id = pp.second->id();

						Vector3d tmp_g(p(1), p(2), p(3));
						double l = proj_normals[face_id].dot(tmp_g);
						Vector3d final_g = tmp_g - l * proj_normals[face_id];
						p(1) = final_g(0);
						p(2) = final_g(1);
						p(3) = final_g(2);
					}*/

					
#endif
					
					

					lnsrch(x, f, g, p, xpls, fpls, mxtake, iretcd, sx, normal_x, normal_y, normal_z);
					if (iretcd == -7) // BAH
					{
						mLastMessage = iretcd;
						return -1;
					}
				}
				else
				{
					dlt = dltsav;
					if (method == 2)
					{
						// SOLVE FOR NEWTON STEP:  AP=-G
						for (i = 1; i <= n; i++)
						{
							wrk1(i) = -g(i);
						}

						lltslv(a, p, wrk1);

						dogdrv(x, f, g, a, p, xpls, fpls, sx, dlt, iretcd, mxtake, wrk0, wrk1, wrk2, wrk3);
					}
					else
					{
						amu = amusav;
						dltp = dlpsav;
						phi = phisav;
						phip0 = phpsav;

						chlhsn(a, sx, udiag);

						// SOLVE FOR NEWTON STEP:  AP=-G
						for (i = 1; i <= n; i++)
						{
							wrk1(i) = -g(i);
						}

						lltslv(a, p, wrk1);

						hookdr(x, f, g, a, udiag, p, xpls, fpls, sx, dlt, iretcd, mxtake, amu, dltp, phi,
							phip0, wrk0, wrk1, wrk2, itncnt);
					}
				}
			}
			// CALCULATE STEP FOR OUTPUT
			for (i = 1; i <= n; i++)
			{
				p(i) = xpls(i) - x(i);
			}

			// CALCULATE GRADIENT AT XPLS
			if (iagflg == -1)
			{
				// CENTRAL DIFFERENCE GRADIENT
				fstocd(xpls, sx, rnf, gpls);
			}
			else if (iagflg == 0)
			{
				// FORWARD DIFFERENCE GRADIENT
				//fstofd(xpls, fpls, gpls, sx, rnf);
				fstofd(xpls, fpls, gpls, sx, rnf, normal_x, normal_y, normal_z);
			}
			else
			{
				// ANALYTIC GRADIENT
				minclass->gradient(xpls, gpls);
			}

			

			// CHECK WHETHER STOPPING CRITERIA SATISFIED
			optstp(xpls, fpls, gpls, x, itncnt, icscmx, itrmcd, sx, fscale, itnlim, iretcd, mxtake, msg);

			if (itrmcd == 0)
			{
				// EVALUATE HESSIAN AT XPLS
				if (iexp != 0)
				{
					if (method == 3)
					{
						secunf(x, g, a, udiag, xpls, gpls, itncnt, rnf, iagflg, noupdt, wrk1, wrk2, wrk3);
					}
					else
					{
						secfac(x, g, a, xpls, gpls, itncnt, rnf, iagflg, noupdt, wrk0, wrk1, wrk2, wrk3);
					}
				}
				else
				{
					if (iahflg == 1)
					{
						minclass->hessian(xpls, a);
					}
					else
					{
						if (iagflg == 1)
						{
							fstofd(xpls, gpls, a, sx, rnf, wrk1);
						}
						else
						{
							//sndofd(xpls, fpls, a, sx, rnf, wrk1, wrk2);
							sndofd(xpls, fpls, a, sx, rnf, wrk1, wrk2, normal_x, normal_y, normal_z);
						}
					}
				}

				if (msg & 8)
				{
					result(xpls, fpls, gpls, a, p, itncnt, 1);
				}

				// X <-- XPLS  AND  G <-- GPLS  AND  F <-- FPLS
				f = fpls;


				//need to be changed here


#if USING_PROJ
				for (i = 1; i <= n; i++)
				{
					
					//x(i) = xpls(i);
					g(i) = gpls(i);
				}

				//project x on the plane

				if (normal_x > 0)
				{
					Point3 query(xpls(1), xpls(2), xpls(3));

					Point3 pp = search_tree->closest_point(query);
					//face_id = pp.second->id();
					//closet point
					x(1) = pp[0];
					x(2) = pp[1];
					x(3) = pp[2];

					//change xpls here
					
					/*xpls(1) = pp[0];
					xpls(2) = pp[1];
					xpls(3) = pp[2];*/
				}
				else
				{
					for (i = 1; i <= n; i++)
					{

						x(i) = xpls(i);
						//g(i) = gpls(i);
					}
				}
#else
				for (i = 1; i <= n; i++)
				{

					x(i) = xpls(i);
					g(i) = gpls(i);
				}
#endif

				


				
			}
		}

		// TERMINATION
		// -----------
		// RESET XPLS,FPLS,GPLS,  IF PREVIOUS ITERATE SOLUTION
		//
		if (itrmcd == 3)
		{
			//itrmcd =  3:    Lower point cannot be found, X is probably optimal
			fpls = f;

			for (i = 1; i <= n; i++)
			{
				xpls(i) = x(i);
				gpls(i) = g(i);
			}
		}
	}

	// PRINT RESULTS
	if (!(msg & 4))
	{
		result(xpls, fpls, gpls, a, p, itncnt, 0);
	}
	mLastMessage = itrmcd;
	return (itrmcd >= 0 && itrmcd <= 3) ? 0 : 1;
}

/*!
  \brief
  Version of Minimize in which the hessian is allocated locally rather than passed
  as an argument to the function.

  \section template_args Template Parameters
  
  \param M  SCPPNT Matrix type.
  \param V  SCPPNT Vector type.
  \param FUNC Type of function to minimize.
  
  \section function_args Function Parameters
   
  \param[in]  &start  Address of vector containing initial estimate of location of minimum.
  \param[out]  &xpls  Address of vector containing local minimum (size must be n on input).
  \param[out]  &fpls  Address of function value at local minimum.
  \param[out]  &gpls  Address of vector containing gradient at local minimum (size must be 
      n on input). 
 */
template<class V, class M, class FUNC> int Uncmin<V, M, FUNC>::Minimize(V &start, V &xpls,
    double &fpls, V &gpls)
{

  M *a = new M(n,n);

  int result = Minimize(start, xpls, fpls, gpls, *a);

  delete a;

  return result;
}

/*!
  \brief 
  Method to set default values for each input variable to the minimization algorithm.
 	
  Translated by Steve Verrill, August 4, 1998.

  \section template_args Template Parameters
  
  \param M  SCPPNT Matrix type.
  \param V  SCPPNT Vector type.
  \param FUNC Type of function to minimize.
 */
template<class V, class M, class FUNC> void Uncmin<V, M, FUNC>::dfault()
{
  // Set default algorithm used for minimization (line search)
  algorithm = 1;

  // SET TYPICAL SIZE OF X AND MINIMIZATION FUNCTION
  typsiz.newsize(n);
  typsiz = 1.0;

  fscale = 1.0;

  // SET TOLERANCES
  mdlt = -1.0;

  SetTolerances(-1.0, -1.0, -1.0);

  stepmx = 0.0;

  // SET FLAGS and options
  iexp = 0;
  ndigit = -1;
  itnlim = 150;
  iagflg = 0;
  iahflg = 0;

  fCheckGradient = 0;
  fCheckHessian = 0;
  fPrintResults = 0;
  fPrintIterationResults = 0;
}

/*!
  \brief 
  Method to check the input for reasonableness.

  Checks for reasonableness of options requested.
  Returns 0 if no serious errors found, otherwise returns 1.

  Translated by Steve Verrill, May 12, 1998.
  
  \section template_args Template Parameters
  
  \param M  SCPPNT Matrix type.
  \param V  SCPPNT Vector type.
  \param FUNC Type of function to minimize.
  
  \section function_args Function Parameters
   
  \param[in]      &x      Address of vector containing initial estimate of location of minimum.
  \param[out]     &sx     Address of vector containing scaling factors.
  \param[in,out]  &msg    Address of error code: If function return value is 1, msg indicates the error found.
  \param[in]      &method Address of message or error code: On input, if the method key includes a term 8, then output will be suppressed.
      On exit, if function return value = 1, then method contains the key of the error message.
 */
template<class V, class M, class FUNC> int Uncmin<V, M, FUNC>::optchk(V &x, V &sx, int &msg,
    int &method)
{
  int i;
  double stpsiz;

  // CHECK THAT PARAMETERS ONLY TAKE ON ACCEPTABLE VALUES.
  // IF NOT, SET THEM TO DEFAULT VALUES.
  bool printMessages = !(msg & 4);

  if (method< 1 || method> 3)method = 1;
  if (iagflg != 1 && iagflg != 0) iagflg = 1;
  if (iahflg != 1 && iahflg != 0) iahflg = 1;
  if (iexp != 0 && iexp != 1) iexp = 1;

  if (!(msg & 1) && iagflg == 0)
  {
    if (printMessages)
    {
      std::fprintf(mfile, "\n\nOPTCHK   User requests that analytic gradient");
      std::fprintf(mfile, " be accepted as properly coded,\n");
      std::fprintf(mfile, "OPTCHK   but an analytic gradient is not");
      std::fprintf(mfile, " supplied,\n");
      std::fprintf(mfile, "OPTCHK   msg = %d,\n", msg);
      std::fprintf(mfile, "OPTCHK   iagflg = %d.\n\n", iagflg);
    }
    msg = -1;
    return 1;
  }

  if (!(msg & 2) && iahflg == 0)
  {
    if (printMessages)
    {
      std::fprintf(mfile, "\n\nOPTCHK   User requests that analytic Hessian");
      std::fprintf(mfile, " be accepted as properly coded,\n");
      std::fprintf(mfile, "OPTCHK   but an analytic Hessian is not");
      std::fprintf(mfile, " supplied,\n");
      std::fprintf(mfile, "OPTCHK   msg = %d,\n", msg);
      std::fprintf(mfile, "OPTCHK   iahflg = %d.\n\n", iahflg);
    }
    msg = -2;
    return 1;
  }

  // CHECK DIMENSION OF PROBLEM

  if (n <= 0)
  {
    if (printMessages) std::fprintf(mfile, "\n\nOPTCHK   Illegal dimension, n = %d\n\n", n);

    msg = -3;
    return 1;
  }

  if (n == 1 && printMessages)
  {
    std::fprintf(mfile, "\n\nOPTCHK   !!!WARNING!!!  This class is ");
    std::fprintf(mfile, "inefficient for problems of size 1.\n");
    std::fprintf(mfile, "OPTCHK   You might want to use a more appropriate routine.\n\n");
  }

  // COMPUTE SCALE MATRIX
  for (i = 1; i <= n; i++)
  {
    if (typsiz(i) == 0) typsiz(i) = 1.0;
    if (typsiz(i) < 0.0) typsiz(i) = -typsiz(i);
    sx(i) = 1.0/typsiz(i);
  }

  // CHECK MAXIMUM STEP SIZE
  if (stepmx <= 0.0)
  {
    stpsiz = 0.0;

    for (i = 1; i <= n; i++)
    {
      stpsiz += x(i)*x(i)*sx(i)*sx(i);
    }

    stpsiz = std::sqrt(stpsiz);

    stepmx = max(1000.0*stpsiz,1000.0);
  }

  // CHECK FUNCTION SCALE
  if (fscale == 0) fscale = 1.0;
  if (fscale < 0.0) fscale = -fscale;

  // CHECK GRADIENT TOLERANCE
  if (gradtl < 0.0)
  {
    if (printMessages) std::fprintf(mfile, "\n\nOPTCHK   Illegal tolerance, gradtl = %lf\n\n", gradtl);

    msg = -4;
    return 1;
  }

  // CHECK ITERATION LIMIT
  if (itnlim < 0)
  {
    if (printMessages)
    {
      std::fprintf(mfile, "\n\nOPTCHK   Illegal iteration limit,");
      std::fprintf(mfile, " itnlim = %d\n\n", itnlim);
    }
    msg = -5;
    return 1;
  }

  // CHECK NUMBER OF DIGITS OF ACCURACY IN FUNCTION FCN
  if (ndigit == 0)
  {
    if (printMessages)
    {
      std::fprintf(mfile, "\n\nOPTCHK   Minimization function has no good");
      std::fprintf(mfile, " digits, ndigit = %d\n\n", ndigit);
    }

    msg = -6;
    return 1;
  }

  if (ndigit < 0) ndigit = (int)(-std::log(epsm)/std::log(10.0));

  // CHECK TRUST REGION RADIUS
  if (mdlt <= 0.0) mdlt = -1.0;
  if (mdlt > stepmx) mdlt = stepmx;

  return 0;
}

/*!
  \brief
  Calculate the first order finite difference approximations for gradients.
  
  Translated by Steve Verrill, April 22, 1998.

  \section template_args Template Parameters
  
  \param M  SCPPNT Matrix type.
  \param V  SCPPNT Vector type.
  \param FUNC Type of function to minimize.
  
  \section function_args Function Parameters
   
  \param[in]      &xpls   Address of new iterate (X[k]).
  \param[in]      fpls    Value of the function to minimize at the new iterate.
  \param[in,out]  &g      Address to vector containing finite difference approximation to the gradient.
  \param[in]      &sx     Address of scaling vector for x.
  \param[in]      rnoise  Relative noise in the function to be minimized.
 	*/                                                  
template<class V, class M, class FUNC> void Uncmin<V, M, FUNC>::fstofd(V &xpls, double fpls, V &g,
    V &sx, double rnoise)                            
{                                                    
  double xmult, stepsz, xtmpj, fhat;                 
  int j;                                             
                                                     
  xmult = std::sqrt(rnoise);

  // gradient

  for (j = 1; j <= n; j++)
  {
    stepsz = xmult*max(std::fabs(xpls(j)), 1.0/sx(j));
    xtmpj = xpls(j);
    xpls(j) = xtmpj + stepsz;

    fhat = minclass->f_to_minimize(xpls);

    xpls(j) = xtmpj;

    g(j) = (fhat - fpls)/stepsz;
  }
}


template<class V, class M, class FUNC> void Uncmin<V, M, FUNC>::fstofd(V &xpls, double fpls, V &g,
	V &sx, double rnoise, int normal_x, int normal_y, int normal_z)
{
	double xmult, stepsz, xtmpj, fhat;
	int j;

	xmult = std::sqrt(rnoise);

	// gradient

	for (j = 1; j <= n; j++)
	{
		stepsz = xmult * max(std::fabs(xpls(j)), 1.0 / sx(j));
		xtmpj = xpls(j);
		xpls(j) = xtmpj + stepsz;

		fhat = minclass->f_to_minimize(xpls, normal_x, normal_y, normal_z);

		xpls(j) = xtmpj;

		g(j) = (fhat - fpls) / stepsz;
	}
}





/*!
  \brief
  Implementation of the fstofd method which finds a finite difference approximation to the Hessian.
  
  Translated by Steve Verrill, April 22, 1998.

  \section template_args Template Parameters
  
  \param M  SCPPNT Matrix type.
  \param V  SCPPNT Vector type.
  \param FUNC Type of function to minimize.
  
  \section function_args Function Parameters
   
  \param[in,out]  &xpls   Address of new iterate.
  \param[in]      &fpls   Address of the gradient of the function to minimize.
  \param[in,out]  &a      Address to Matrix containing finite difference 
      approximation to the hessian. Only the elements in the diagonal and 
      lower-half are returned.
  \param[in]      &sx     Address of scaling vector for x.
  \param[in]      rnoise  Relative noise in the function to be minimized.
  \param[in]      &fhat   Address of workspace vector.
 */                                                  
template<class V, class M, class FUNC> void Uncmin<V, M, FUNC>::fstofd(V &xpls, V &fpls, M &a,
    V &sx, double rnoise, V &fhat)
{
  double xmult, stepsz, xtmpj;
  int i, j, nm1;

  xmult = std::sqrt(rnoise);

  for (j = 1; j <= n; j++)
  {
    stepsz = xmult * max(std::fabs(xpls(j)), 1.0/sx(j));
    xtmpj = xpls(j);
    xpls(j) = xtmpj + stepsz;

    minclass->gradient(xpls, fhat);

    xpls(j) = xtmpj;

    for (i = 1; i <= n; i++)
    {
      a(i, j) = (fhat(i) - fpls(i))/stepsz;
    }
  }

  nm1 = n - 1;

  for (j = 1; j <= nm1; j++)
  {
    for (i = j+1; i <= n; i++)
    {
      a(i, j) = (a(i,j) + a(j,i))/2.0;
    }
  }
  return;
}

/*!
	\brief 
	
	Method to check the analytic gradient supplied by the user.
	
	This function compares the analytic gradient to the first-order finite difference gradient.
  The function return value zero if gradient checks, or 1 if the analytic gradient appears to be off.
	
  Translated by Steve Verrill, April 22, 1998.
  
  \section template_args Template Parameters

  \param M  SCPPNT Matrix type.
  \param V  SCPPNT Vector type.
  \param FUNC Type of function to minimize.

  \section function_args Function Parameters
 
  \param[in]  &x      Address of location vector at which the gradient is to be checked.
  \param[in]  f       Function value.
  \param[in]  &g      Address of vector containing analytic gradient.
  \param[in]  &typsiz Address of scale vector for x.
  \param[in]  &sx     Relative noise in the function to be minimized.
  \param[in]  fscale  Estimate of scale of f_to_minimize. 
  \param[in]  rnf     Relative noise in f_to_minimize.
  \param[in]  analtl  Tolerance for comparison of estimated and analytical gradient.
  \param[in]  &gest   Address of finite difference gradient
  \param[out] &msg    Address of message or error code: On input, if msg code includes a term 8, suppresses output.
      On output, if function return value = 1, then contains the key of the error message.
 */                                                  
template<class V, class M, class FUNC> int Uncmin<V, M, FUNC>::grdchk(V &x, double f, V &g,
    V &typsiz, V &sx, double fscale, double rnf, double analtl, V &gest, int &msg)
{
  double gs;
  int ker, i;

  // COMPUTE FIRST ORDER FINITE DIFFERENCE GRADIENT AND COMPARE TO
  // ANALYTIC GRADIENT.
  fstofd(x, f, gest, sx, rnf);

  ker = 0;

  for (i = 1; i <= n; i++)
  {
    gs = max(std::fabs(f), fscale)/max(std::fabs(x(i)), typsiz(i));

    if (std::fabs(g(i) - gest(i)) > max(std::fabs(g(i)), gs)*analtl)
      ker = 1;
  }

  if (ker == 0)
    return 0;

  if (!(msg & 4))
  {
    std::fprintf(mfile, "\nThere appears to be an error in the coding");
    std::fprintf(mfile, " of the gradient method.\n\n\n");
    std::fprintf(mfile, "Component   Analytic   Finite Difference\n\n");

    for (i = 1; i <= n; i++)
    {
      std::fprintf(mfile, "%d  %lf  %lf\n", i, (double) g(i), (double) gest(i));
    }
  }
  msg = -21;
  return 1;
}

/*!
  \brief
  Test for convergence of iterations.
  
  Determines whether the algorithm should terminate due to any of the following:
   1. Problem solved within user tolerance; 
   2. Convergence within user tolerance; 
   3. Iteration limit reached; 
   4. Divergence or too restrictive maximum step (stepmx) suspected. 
   
   Translated by Steve Verrill, May 12, 1998.

  \section template_args Template Parameters

  \param M  SCPPNT Matrix type.
  \param V  SCPPNT Vector type.
  \param FUNC Type of function to minimize.

  \section function_args Function Parameters
 
  \param[in]  &xpls   Address of location vector at current iteration (X[k]).
  \param[in]  fpls    Function value at xpls.
  \param[in]  &gpls   Address of gradient or approximation.
  \param[in]  &x      Address of location vector at previous iteration (X[k-1]).
  \param[in]  &itncnt Address of current iteration number.
  \param[out] &icscmx Address of number of consecutive steps >= stepmx (retain 
      between successive calls).
  \param[out] &itrmcd Address of termination code.
  \param[in]  &sx     Address of scaling vector for x.
  \param[in]  &fscale Address of estimate of scale of f_to_minimize. 
  \param[in]  &itnlim Address of maximum number of allowable iterations.
  \param[in]  &iretcd Address of return code.
  \param[in]  &mxtake Address of boolean flag indicating step of maximum length was 
      used.
  \param[out] &msg    Address of message or error code: If msg code includes a term 8, suppresses output.
 */
template<class V, class M, class FUNC> void Uncmin<V, M, FUNC>::optstp(V &xpls, double fpls,
    V &gpls, V &x, int &itncnt, int &icscmx, int &itrmcd, V &sx, double &fscale, int &itnlim,
    int &iretcd, bool &mxtake, int &msg)
{
  int i;
  double d, rgx, relgrd, rsx = 0.0, relstp; // Initialized rsx = 0.0. ww, 12-20-2007

  itrmcd = 0;

  // LAST GLOBAL STEP FAILED TO LOCATE A POINT LOWER THAN X
  if (iretcd == 1)
  {
    itrmcd = 3;

    if (!(msg & 4))
    {
      std::fprintf(mfile, "\n\nOPTSTP    The last global step failed");
      std::fprintf(mfile, " to locate a point lower than x.\n");
      std::fprintf(mfile, "OPTSTP    Either x is an approximate local");
      std::fprintf(mfile, " minimum of the function,\n");
      std::fprintf(mfile, "OPTSTP    the function is too nonlinear for");
      std::fprintf(mfile, " this algorithm, or\n");
      std::fprintf(mfile, "OPTSTP    steptl is too large.\n");
    }
    return;
  }
  else
  {
    // FIND DIRECTION IN WHICH RELATIVE GRADIENT MAXIMUM.
    d = max(std::fabs(fpls), fscale);
    rgx = 0.0;

    for (i = 1; i <= n; i++)
    {
      relgrd = std::fabs(gpls(i))*max(std::fabs(xpls(i)), 1.0/sx(i))/d;
      rgx = max(rgx, relgrd);
    }
    mLastGradCrit = rgx; // BAH

    // FIND DIRECTION IN WHICH RELATIVE STEPSIZE MAXIMUM
    // Moved before check of gradient criterion so it is always computed - BAH
    if (itncnt > 0)
    {
      rsx = 0.0;

      for (i = 1; i <= n; i++)
      {
        relstp = std::fabs(xpls(i) - x(i))/ max(std::fabs(xpls(i)), 1.0/sx(i));

        rsx = max(rsx, relstp);
      }
      mLastStepCrit = rsx;
    }
    // Check whether gradient criterion is met
    if (rgx <= gradtl)
    {
      itrmcd = 1;

      if (!(msg & 4))
      {
        std::fprintf(mfile, "\n\nOPTSTP    The relative gradient is close");
        std::fprintf(mfile, " to zero.\n");
        std::fprintf(mfile, "OPTSTP    The current iterate is probably");
        std::fprintf(mfile, " a solution.\n");
      }
      return;
    }
    if (itncnt == 0)
      return;

    // Check whether maximum step size criteria is met
    if (rsx <= steptl)
    {
      itrmcd = 2;

      if (!(msg & 4))
      {
        std::fprintf(mfile, "\n\nOPTSTP    Successive iterates are within");
        std::fprintf(mfile, " steptl.\n");
        std::fprintf(mfile, "OPTSTP    The current iterate is probably");
        std::fprintf(mfile, " a solution.\n");
      }
      return;
    }

    // CHECK ITERATION LIMIT
    if (itncnt >= itnlim)
    {
      itrmcd = 4;

      if (!(msg & 4))
      {
        std::fprintf(mfile, "\n\nOPTSTP    The iteration limit was reached.\n");
        std::fprintf(mfile, "OPTSTP    The algorithm failed.\n");
      }

      return;
    }

    // CHECK NUMBER OF CONSECUTIVE STEPS \ STEPMX
    if (!mxtake)
    {
      icscmx = 0;

      return;
    }

    if (!(msg & 4))
    {
      std::fprintf(mfile, "\n\nOPTSTP    Step of maximum length (stepmx)");
      std::fprintf(mfile, " taken.\n");
    }

    icscmx++;

    if (icscmx < 5)
      return;

    itrmcd = 5;

    if (!(msg & 4))
    {
      std::fprintf(mfile, "\n\nOPTSTP    Maximum step size exceeded");
      std::fprintf(mfile, " five consecutive times.\n");
      std::fprintf(mfile, "OPTSTP    Either the function is unbounded");
      std::fprintf(mfile, " below,\n");
      std::fprintf(mfile, "OPTSTP    becomes asymptotic to a finite value");
      std::fprintf(mfile, " from above in some direction, or\n");
      std::fprintf(mfile, "OPTSTP    stepmx is too small.\n");
    }
  }
}

/*!
  \brief
  Provides the initial Hessian when secant updates are being used.
  
  Translated by Steve Verrill, April 27, 1998.

  \section template_args Template Parameters

  \param M  SCPPNT Matrix type.
  \param V  SCPPNT Vector type.
  \param FUNC Type of function to minimize.

  \section function_args Function Parameters
 
  \param[out] &a      Address of initial hessian (lower triangular matrix)
  \param[in]  &sx     Address of scaling vector for x.
  \param[in]  method  Indictor for algorithm to use to solve the minimization problem.
               Method = 1,2: Use factored secant method;
               method = 3:  Use unfactored secant method.
 */

template<class V, class M, class FUNC> void Uncmin<V, M, FUNC>::hsnint(M &a, V &sx, int method)
{
  int i, j;

  for (j = 1; j <= n; j++)
  {
    if (method == 3)
    {
      a(j, j) = sx(j)*sx(j);
    }
    else
    {
      a(j, j) = sx(j);
    }

    for (i = j + 1; i <= n; i++)
    {
      a(i, j) = 0.0;
    }
  }
  return;
}

/*!
  \brief
  Finds the second order forward finite difference approximations to the Hessian.  
  
  For optimization use this method to estimate the Hessian of the optimization function
  if no analytical user function has been supplied for either the gradient or the Hessian, 
  and the optimization function is inexpensive to evaluate.
  
  Translated by Steve Verrill, May 8, 1998.

  \section template_args Template Parameters

  \param M  SCPPNT Matrix type.
  \param V  SCPPNT Vector type.
  \param FUNC Type of function to minimize.

  \section function_args Function Parameters
 
  \param[out]   &xpls       Address of new iterate, (X[k]).
  \param[in]    &fpls       Address of function value at the new iterate.
  \param[out]   &a          Address of finite difference approximation to the hessian. Only the entries in
      the lower triangular matrix and diagonal are returned.
  \param[in]    &sx         Address of scaling vector for x.
  \param[in]    &rnoise     Address of relative noise in the function to be minimized.
  \param[out]   &stepsz     Address of workspace (stepsize in i-th component direction).
  \param[out]   &anbr       Address of workspace (neighbor in i-th direction).
 */
template<class V, class M, class FUNC> void Uncmin<V, M, FUNC>::sndofd(V &xpls, double &fpls, M &a,
    V &sx, double &rnoise, V &stepsz, V &anbr)
{
  double xmult, xtmpi, xtmpj, fhat;
  int i, j;

  // FIND I-TH STEPSIZE AND EVALUATE NEIGHBOR IN DIRECTION
  // OF I-TH UNIT VECTOR
  xmult = std::pow(rnoise, 1.0/3.0);

  for (i = 1; i <= n; i++)
  {
    stepsz(i) = xmult*max(std::fabs(xpls(i)),1.0/sx(i));
    xtmpi = xpls(i);
    xpls(i) = xtmpi + stepsz(i);

    anbr(i) = minclass->f_to_minimize(xpls);
    xpls(i) = xtmpi;
  }

  // CALCULATE COLUMN I OF A
  for (i = 1; i <= n; i++)
  {
    xtmpi = xpls(i);
    xpls(i) = xtmpi + 2.0*stepsz(i);
    fhat = minclass->f_to_minimize(xpls);
    a(i, i) = ((fpls - anbr(i)) + (fhat - anbr(i)))/
    (stepsz(i)*stepsz(i));

    // CALCULATE SUB-DIAGONAL ELEMENTS OF COLUMN
    if (i != n)
    {
      xpls(i) = xtmpi + stepsz(i);

      for (j = i+1; j <= n; j++)
      {
        xtmpj = xpls(j);
        xpls(j) = xtmpj + stepsz(j);
        fhat = minclass->f_to_minimize(xpls);
        a(j, i) = ((fpls - anbr(i)) + (fhat - anbr(j)))/
        (stepsz(i)*stepsz(j));
        xpls(j) = xtmpj;
      }
    }
    xpls(i) = xtmpi;
  }
}


template<class V, class M, class FUNC> void Uncmin<V, M, FUNC>::sndofd(V &xpls, double &fpls, M &a,
	V &sx, double &rnoise, V &stepsz, V &anbr, int normal_x, int normal_y, int normal_z)
{
	double xmult, xtmpi, xtmpj, fhat;
	int i, j;

	// FIND I-TH STEPSIZE AND EVALUATE NEIGHBOR IN DIRECTION
	// OF I-TH UNIT VECTOR
	xmult = std::pow(rnoise, 1.0 / 3.0);

	for (i = 1; i <= n; i++)
	{
		stepsz(i) = xmult * max(std::fabs(xpls(i)), 1.0 / sx(i));
		xtmpi = xpls(i);
		xpls(i) = xtmpi + stepsz(i);

		//anbr(i) = minclass->f_to_minimize(xpls);
		anbr(i) = minclass->f_to_minimize(xpls, normal_x, normal_y, normal_z);
		xpls(i) = xtmpi;
	}

	// CALCULATE COLUMN I OF A
	for (i = 1; i <= n; i++)
	{
		xtmpi = xpls(i);
		xpls(i) = xtmpi + 2.0*stepsz(i);
		//fhat = minclass->f_to_minimize(xpls);
		fhat = minclass->f_to_minimize(xpls, normal_x, normal_y, normal_z);
		a(i, i) = ((fpls - anbr(i)) + (fhat - anbr(i))) /
			(stepsz(i)*stepsz(i));

		// CALCULATE SUB-DIAGONAL ELEMENTS OF COLUMN
		if (i != n)
		{
			xpls(i) = xtmpi + stepsz(i);

			for (j = i + 1; j <= n; j++)
			{
				xtmpj = xpls(j);
				xpls(j) = xtmpj + stepsz(j);
				//fhat = minclass->f_to_minimize(xpls);
				fhat = minclass->f_to_minimize(xpls, normal_x, normal_y, normal_z);
				a(j, i) = ((fpls - anbr(i)) + (fhat - anbr(j))) /
					(stepsz(i)*stepsz(j));
				xpls(j) = xtmpj;
			}
		}
		xpls(i) = xtmpi;
	}
}


/*!
  \brief
  Checks the analytic hessian supplied by the user.
  
  Returns 0 if hessian checks, otherwise returns 1.
  
  Translated by Steve Verrill, April 23, 1998.

  \section template_args Template Parameters

  \param M  SCPPNT Matrix type.
  \param V  SCPPNT Vector type.
  \param FUNC Type of function to minimize.

  \section function_args Function Parameters

  \param[in]  &x          Address of iterate (X[k]) at which the Hessian is to be checked.
  \param[in]  &f          Address of function value f(x).
  \param[out]  &g         Address of gradient g(x).
  \param[out]  &a         Address of hessian h(x). On exit: Hessian in lower triangle.
  \param[in]  &typsiz     Address of typical size for each component of x.
  \param[in]  &sx         Address of scaling vector for x:  sx[i] = 1.0/typsiz[i].
  \param[in]  rnf         Relative noise in f_to_minimize.
  \param[in]  analtl      Tolerance for comparison of estimated and analytic gradients.
  \param[in]  &iagflg     Address of flag: If iagflg = 1, then an analytic gradient is being supplied.
  \param[in]  &udiag      Address of workspace vector.
  \param[in]  &wrk1       Address of workspace vector.
  \param[in]  &wrk2       Address of workspace vector.
  \param[in,out]  &msg    Address of message or error code. On input, if msg code contains factor 8, then print 
      discrepancy report for analytic hessian. On output, msg = -22 indicates probably coding error
      of hessian.
 */
template<class V, class M, class FUNC> int Uncmin<V, M, FUNC>::heschk(V &x, double &f, V &g, M &a,
    V &typsiz, V &sx, double rnf, double analtl, int &iagflg, V &udiag, V &wrk1, V &wrk2, int &msg)
{
  int i, j, ker;
  double hs;

  // COMPUTE FINITE DIFFERENCE APPROXIMATION H TO THE HESSIAN.
  if (iagflg == 1)
    fstofd(x, g, a, sx, rnf, wrk1);
  else
    sndofd(x, f, a, sx, rnf, wrk1, wrk2);

  ker = 0;

  // COPY LOWER TRIANGULAR PART OF H TO UPPER TRIANGULAR PART
  // AND DIAGONAL OF H TO UDIAG
  for (j = 1; j <= n; j++)
  {
    udiag(j) = a(j,j);

    for (i = j + 1; i <= n; i++)
    {
      a(j, i) = a(i,j);
    }
  }

  // COMPUTE ANALYTIC HESSIAN AND COMPARE TO FINITE DIFFERENCE
  // APPROXIMATION.
  minclass->hessian(x, a);

  for (j = 1; j <= n; j++)
  {
    hs = max(std::fabs(g(j)), 1.0)/max(std::fabs(x(j)), typsiz(j));

    if (std::fabs(a(j, j) - udiag(j)) > max(std::fabs(udiag(j)), hs)*analtl)
      ker = 1;

    for (i = j + 1; i <= n; i++)
    {
      if (std::fabs(a(i, j) - a(j, i)) > max(std::fabs(a(i, j)), hs)*analtl)
        ker = 1;
    }
  }
  if (ker == 0)
    return 0;

  if (!(msg & 4))
  {
    std::fprintf(mfile, "\nThere appears to be an error in the coding");
    std::fprintf(mfile, " of the Hessian method.\n\n\n");
    std::fprintf(mfile, "Row   Column   Analytic   Finite Difference\n\n");

    for (i = 1; i <= n; i++)
    {
      for (j = 1; j < i; j++)
      {
        std::fprintf(mfile, "%d  %d  %lf  %lf\n", i, j, (double) a(i, j), (double) a(j, i));
      }
      std::fprintf(mfile, "%d  %d  %lf  %lf\n", i, i, (double) a(i, i), (double) udiag(i));
    }
  }
  msg = -22;
  return 1;
}


template<class V, class M, class FUNC> int Uncmin<V, M, FUNC>::heschk(V &x, double &f, V &g, M &a,
	V &typsiz, V &sx, double rnf, double analtl, int &iagflg, V &udiag, V &wrk1, V &wrk2, int &msg, int normal_x, int normal_y, int normal_z)
{
	int i, j, ker;
	double hs;

	// COMPUTE FINITE DIFFERENCE APPROXIMATION H TO THE HESSIAN.
	if (iagflg == 1)
		fstofd(x, g, a, sx, rnf, wrk1);
	else
		sndofd(x, f, a, sx, rnf, wrk1, wrk2, normal_x, normal_y, normal_z);

	ker = 0;

	// COPY LOWER TRIANGULAR PART OF H TO UPPER TRIANGULAR PART
	// AND DIAGONAL OF H TO UDIAG
	for (j = 1; j <= n; j++)
	{
		udiag(j) = a(j, j);

		for (i = j + 1; i <= n; i++)
		{
			a(j, i) = a(i, j);
		}
	}

	// COMPUTE ANALYTIC HESSIAN AND COMPARE TO FINITE DIFFERENCE
	// APPROXIMATION.
	minclass->hessian(x, a);

	for (j = 1; j <= n; j++)
	{
		hs = max(std::fabs(g(j)), 1.0) / max(std::fabs(x(j)), typsiz(j));

		if (std::fabs(a(j, j) - udiag(j)) > max(std::fabs(udiag(j)), hs)*analtl)
			ker = 1;

		for (i = j + 1; i <= n; i++)
		{
			if (std::fabs(a(i, j) - a(j, i)) > max(std::fabs(a(i, j)), hs)*analtl)
				ker = 1;
		}
	}
	if (ker == 0)
		return 0;

	if (!(msg & 4))
	{
		std::fprintf(mfile, "\nThere appears to be an error in the coding");
		std::fprintf(mfile, " of the Hessian method.\n\n\n");
		std::fprintf(mfile, "Row   Column   Analytic   Finite Difference\n\n");

		for (i = 1; i <= n; i++)
		{
			for (j = 1; j < i; j++)
			{
				std::fprintf(mfile, "%d  %d  %lf  %lf\n", i, j, (double)a(i, j), (double)a(j, i));
			}
			std::fprintf(mfile, "%d  %d  %lf  %lf\n", i, i, (double)a(i, i), (double)udiag(i));
		}
	}
	msg = -22;
	return 1;
}

/*!
  \brief
  Prints information.
  
  Translated by Steve Verrill, May 11, 1998.

  \section template_args Template Parameters

  \param M  SCPPNT Matrix type.
  \param V  SCPPNT Vector type.
  \param FUNC Type of function to minimize.

  \section function_args Function Parameters

  \param[in]  &x       Address of iterate (X[k]) at which the Hessian is to be checked.
  \param[in]  &f       Address of function value f(x).
  \param[in]  &g       Address of gradient g(x).
  \param[in]  &a       Address of hessian h(x). On exit: Hessian in lower triangle.
  \param[in]  &p       Address of step taken.
  \param[in]  &itncnt  Address of iteration number (k).
  \param[in]  iflg     Flag controlling the information to print
 */
template<class V, class M, class FUNC> void Uncmin<V, M, FUNC>::result(V &x, double &f, V &g, M &a,
    V &p, int &itncnt, int iflg)
{
  int i, j, iii, num5, remain, iii5, iiir;
  int ilow, ihigh;

  num5 = n/5;
  remain = n%5;

  // PRINT ITERATION NUMBER
  std::fprintf(mfile, "\n\nRESULT      Iterate k = %d\n", itncnt);

  if (iflg != 0)
  {
    // PRINT STEP
    std::fprintf(mfile, "\n\nRESULT      Step\n\n");

    ilow = -4;
    ihigh = 0;

    for (i = 1; i <= num5; i++)
    {
      ilow += 5;
      ihigh += 5;

      std::fprintf(mfile, "%d--%d     ", ilow, ihigh);

      for (j = 1; j <= 5; j++)
      {
        std::fprintf(mfile, "%lf  ", (double) p(ilow+j-1));
      }

      std::fprintf(mfile, "\n");
    }

    ilow += 5;
    ihigh = ilow + remain - 1;

    std::fprintf(mfile, "%d--%d     ", ilow, ihigh);

    for (j = 1; j <= remain; j++)
    {
      std::fprintf(mfile, "%lf  ", (double) p(ilow+j-1));
    }

    std::fprintf(mfile, "\n");
  }

  // PRINT CURRENT ITERATE
  std::fprintf(mfile, "\n\nRESULT      Current x\n\n");

  ilow = -4;
  ihigh = 0;

  for (i = 1; i <= num5; i++)
  {
    ilow += 5;
    ihigh += 5;

    std::fprintf(mfile, "%d--%d     ", ilow, ihigh);

    for (j = 1; j <= 5; j++)
    {
      std::fprintf(mfile, "%lf  ", (double) x(ilow+j-1));
    }
    std::fprintf(mfile, "\n");
  }

  ilow += 5;
  ihigh = ilow + remain - 1;

  std::fprintf(mfile, "%d--%d     ", ilow, ihigh);

  for (j = 1; j <= remain; j++)
  {
    std::fprintf(mfile, "%lf  ", (double) x(ilow+j-1));
  }

  std::fprintf(mfile, "\n");

  // PRINT FUNCTION VALUE
  std::fprintf(mfile, "\n\nRESULT      f_to_minimize at x = %lf\n", f);

  // PRINT GRADIENT
  std::fprintf(mfile, "\n\nRESULT      Gradient at x\n\n");

  ilow = -4;
  ihigh = 0;

  for (i = 1; i <= num5; i++)
  {
    ilow += 5;
    ihigh += 5;

    std::fprintf(mfile, "%d--%d     ", ilow, ihigh);

    for (j = 1; j <= 5; j++)
    {
      std::fprintf(mfile, "%lf  ", (double) g(ilow+j-1));
    }
    std::fprintf(mfile, "\n");
  }

  ilow += 5;
  ihigh = ilow + remain - 1;

  std::fprintf(mfile, "%d--%d     ", ilow, ihigh);

  for (j = 1; j <= remain; j++)
  {
    std::fprintf(mfile, "%lf  ", (double) g(ilow+j-1));
  }

  std::fprintf(mfile, "\n");

  // PRINT HESSIAN FROM ITERATION K
  if (iflg != 0)
  {
    std::fprintf(mfile, "\n\nRESULT      Hessian at x\n\n");

    for (iii = 1; iii <= n; iii++)
    {
      iii5 = iii/5;
      iiir = iii%5;

      ilow = -4;
      ihigh = 0;

      for (i = 1; i <= iii5; i++)
      {
        ilow += 5;
        ihigh += 5;

        std::fprintf(mfile, "i = %d, j = ", iii);
        std::fprintf(mfile, "%d--%d     ", ilow, ihigh);

        for (j = 1; j <= 5; j++)
        {
          std::fprintf(mfile, "%lf  ", (double) a(iii, ilow+j-1));
        }
        std::fprintf(mfile, "\n");
      }
      ilow += 5;
      ihigh = ilow + iiir - 1;

      std::fprintf(mfile, "i = %d, j = ", iii);
      std::fprintf(mfile, "%d--%d     ", ilow, ihigh);

      for (j = 1; j <= iiir; j++)
      {
        std::fprintf(mfile, "%lf  ", (double) a(iii, ilow+j-1));
      }

      std::fprintf(mfile, "\n");
    }
  }
}

/*!
  \brief
  Solves Ax = b where A is an upper triangular matrix.
    
  Note that A is input as a lower triangular matrix and this method takes its transpose 
  implicitly.

  Translated by Steve Verrill, April 14, 1998.

  \section template_args Template Parameters

  \param M  SCPPNT Matrix type.
  \param V  SCPPNT Vector type.
  \param FUNC Type of function to minimize.

  \section function_args Function Parameters

  \param[in]  &a         Address of n by n lower triangular matrix A (preserved).
  \param[out] &x         Address of solution vector x.
  \param[in]  &b         Address of right-hand side vector b.
 */
template<class V, class M, class FUNC> void Uncmin<V, M, FUNC>::bakslv(M &a, V &x, V &b)
{
  int i, ip1, j;
  double sum;

  // SOLVE (L-TRANSPOSE)X=B. (BACK SOLVE)
  i = n;
  x(i) = b(i)/a(i,i);

  while (i > 1)
  {
    ip1 = i;
    i--;

    sum = 0.0;

    for (j = ip1; j <= n; j++)
    {
      sum += a(j, i)*x(j);
    }
    x(i) = (b(i) - sum)/a(i,i);
  }
}

/*!
  \brief
  Finds the LL' decomposition of the perturbed model hessian matrix
  A+mu*I (where mu > 0 and I is the identity matrix) which is safely
  positive definite. If A is strictly positiv definite upon entry,
  then mu = 0.
  
  Translated by Steve Verrill, April 14, 1998.

  Description: 
  
  1. If A has any negative diagonal elements, then choose mu > 0 such
  that the diagonal of A1 = A + mu*I is all positive, with the ratio
  of its smallest to largest elements in the order of sqrt(epsm).
  
  2. A1 undergoes a perturbed Cholesky decomposition which results in
  and LL' decomposition of A1 + D, where D is a non-negative diagonal
  matrix which is implicitly added to A1 during the decomposition, if
  A1 is not positive definite. A1 is retained and not changed during
  this process by copying L into the upper triangular part of A1 and
  the diagonal into udiag. Then the Cholesky decomposition routine is
  called. On return, addmax contains the maximum element of D.
  
  3. If addmax = 0, A1 was positive definite going into step 2 and
  the process returns to the calling program. Otherwise, the minimum
  number sdd which must be added to the diagonal of A to make it safely
  strictly diagonallly dominant is calculated. Since A + addmax*I and
  A + sdd*I are both safely positive definite, choose 
  mu = min( addmax, sdd) and decompose A + mu*I to obtain L.

  \section template_args Template Parameters

  \param M  SCPPNT Matrix type.
  \param V  SCPPNT Vector type.
  \param FUNC Type of function to minimize.

  \section function_args Function Parameters

  \param[in,out]  &a      Address of nxn matrix A: On entry, A is the model hessian 
      (only the lower triangle and diagonal stored). 
      On exit, A contains the L factor of the LL' decomposition of the perturbed model hessian 
      in the lower triangle and diagonal, and contains the off-diagonal elements of the 
      hessian in the upper triangle (Note: diagonal elements of the hessian are returned 
      in udiag).
  \param[in]      &sx     Address of scaling vector for x.
  \param[out]     &udiag  Address of diagonal of the hessian (returned on exit).              
 */
template<class V, class M, class FUNC> void Uncmin<V, M, FUNC>::chlhsn(M &a, V &sx, V &udiag)
{
  int i, j, im1, jm1;
  double tol, diagmx, diagmn, posmax, amu, offmax, evmin, evmax, offrow, sdd;

  double addmax;

  // SCALE HESSIAN
  // PRE- AND POST- MULTIPLY "A" BY INV(SX)
  for (j = 1; j <= n; j++)
  {
    for (i = j; i <= n; i++)
    {
      a(i, j) /= (sx(i)*sx(j));
    }
  }

  // STEP1
  // -----
  // NOTE:  IF A DIFFERENT TOLERANCE IS DESIRED THROUGHOUT THIS
  // ALGORITHM, CHANGE TOLERANCE HERE:
  tol = std::sqrt(epsm);

  diagmx = a(1, 1);
  diagmn = a(1, 1);

  for (i = 2; i <= n; i++)
  {
    if (a(i, i) < diagmn)
      diagmn = a(i, i);
    if (a(i, i) > diagmx)
      diagmx = a(i, i);
  }
  posmax = max(diagmx, 0.0);

  // DIAGMN .LE. 0
  if (diagmn <= posmax*tol)
  {
    amu = tol*(posmax - diagmn) - diagmn;

    if (amu == 0.0)
    {
      // FIND LARGEST OFF-DIAGONAL ELEMENT OF A
      offmax = 0.0;

      for (i = 2; i <= n; i++)
      {
        im1 = i - 1;

        for (j = 1; j <= im1; j++)
        {
          if (std::fabs(a(i, j)) > offmax)
            offmax = std::fabs(a(i, j));
        }
      }
      amu = offmax;

      if (amu == 0.0)
      {
        amu = 1.0;
      }
      else
      {
        amu *= 1.0 + tol;
      }
    }

    // A = A + MU*I
    for (i = 1; i <= n; i++)
    {
      a(i, i) += amu;
    }

    diagmx += amu;
  }

  // STEP2
  // -----
  // COPY LOWER TRIANGULAR PART OF "A" TO UPPER TRIANGULAR PART
  // AND DIAGONAL OF "A" TO UDIAG
  for (j = 1; j <= n; j++)
  {
    udiag(j) = a(j,j);
    for (i = j + 1; i <= n; i++)
    {
      a(j, i) = a(i,j);
    }
  }
  choldc(a, diagmx, tol, addmax);

  // STEP3
  // -----
  // IF ADDMAX = 0, "A" WAS POSITIVE DEFINITE GOING INTO STEP 2,
  // THE LL' DECOMPOSITION HAS BEEN DONE, AND WE RETURN.
  // OTHERWISE, ADDMAX > 0.  PERTURB "A" SO THAT IT IS SAFELY
  // DIAGONALLY DOMINANT AND FIND LL' DECOMPOSITION
  if (addmax > 0.0)
  {
    // RESTORE ORIGINAL "A" (LOWER TRIANGULAR PART AND DIAGONAL)
    for (j = 1; j <= n; j++)
    {
      a(j, j) = udiag(j);

      for (i = j + 1; i <= n; i++)
      {
        a(i, j) = a(j,i);
      }
    }

    // FIND SDD SUCH THAT A + SDD*I IS SAFELY POSITIVE DEFINITE
    // NOTE:  EVMIN < 0 SINCE A IS NOT POSITIVE DEFINITE;
    evmin = 0.0;
    evmax = a(1, 1);
    for (i = 1; i <= n; i++)
    {
      offrow = 0.0;
      im1 = i - 1;

      for (j = 1; j <= im1; j++)
      {
        offrow += std::fabs(a(i, j));
      }

      for (j = i + 1; j <= n; j++)
      {
        offrow += std::fabs(a(j, i));
      }

      evmin = min(evmin, a(i, i) - offrow);
      evmax = max(evmax, a(i, i) + offrow);
    }
    sdd = tol*(evmax - evmin) - evmin;

    // PERTURB "A" AND DECOMPOSE AGAIN
    amu = min(sdd, addmax);

    for (i = 1; i <= n; i++)
    {
      a(i, i) += amu;
      udiag(i) = a(i,i);
    }

    // "A" NOW GUARANTEED SAFELY POSITIVE DEFINITE
    choldc(a, 0.0, tol, addmax);
  }

  // UNSCALE HESSIAN AND CHOLESKY DECOMPOSITION MATRIX
  for (j = 1; j <= n; j++)
  {
    for (i = j; i <= n; i++)
    {
      a(i, j) *= sx(i);
    }

    jm1 = j - 1;

    for (i = 1; i <= jm1; i++)
    {
      a(i, j) *= sx(i)*sx(j);
    }

    udiag(j) *= sx(j)*sx(j);
  }
  return;
}

/*!
  \brief
  Finds the perturbed LL' decomposition of A + D, where D is a non-
  negative diagonal matrix added to A if necessary to allow the
  Cholesky decomposition to continue.
  
  The normal Cholesky decomposition is performed. However, if at any
  point the algorithm would attempt to set L(i,i) = SQRT(Temp) with
  Temp < Tol*Diagmx, then L(i,i) is set to SQRT(Tol*Diagmx) instead.
  This is equivalent to adding Tol_Diagmx - Temp to A(i,i).
    
  Translated by Steve Verrill, April 15, 1998.

  \section template_args Template Parameters

  \param M  SCPPNT Matrix type.
  \param V  SCPPNT Vector type.
  \param FUNC Type of function to minimize.

  \section function_args Function Parameters

  \param[in,out]  &a      Address of nxn matrix A: On entry, A is the matrix for which 
      to find the perturbed Cholesky decomposition. On exit, A contains L of the LL' 
      decomposition in lower triangle including the diagonal.
  \param[in]      diagmx  Maximum diagonal element of matrix A.
  \param[in]      tol     Tolerance.
  \param[out]     &addmax Address of maximum amount implicitly added to diagonal 
      of A in forming the Cholesky decomposition of A + D.            
 */ 
template<class V, class M, class FUNC> void Uncmin<V, M, FUNC>::choldc(M &a, double diagmx,
    double tol, double &addmax)
{
  int i, j, jm1, jp1, k;
  double aminl, amnlsq, offmax, sum, temp;

  addmax = 0.0;
  aminl = std::sqrt(diagmx*tol);
  amnlsq = aminl*aminl;

  // FORM COLUMN J OF L
  for (j = 1; j <= n; j++)
  {
    // FIND DIAGONAL ELEMENTS OF L
    sum = 0.0;
    jm1 = j - 1;
    jp1 = j + 1;

    for (k = 1; k <= jm1; k++)
    {
      sum += a(j, k)*a(j, k);
    }
    temp = a(j, j) - sum;

    if (temp >= amnlsq)
    {
      a(j, j) = std::sqrt(temp);
    }
    else
    {
      // FIND MAXIMUM OFF-DIAGONAL ELEMENT IN COLUMN
      offmax = 0.0;

      for (i = jp1; i <= n; i++)
      {
        if (std::fabs(a(i, j)) > offmax)
          offmax = std::fabs(a(i, j));
      }
      if (offmax <= amnlsq)
        offmax = amnlsq;

      // ADD TO DIAGONAL ELEMENT  TO ALLOW CHOLESKY DECOMPOSITION TO CONTINUE
      a(j, j) = std::sqrt(offmax);
      addmax = max(addmax, offmax-temp);
    }

    // FIND I,J ELEMENT OF LOWER TRIANGULAR MATRIX
    for (i = jp1; i <= n; i++)
    {
      sum = 0.0;

      for (k = 1; k <= jm1; k++)
      {
        sum += a(i, k)*a(j, k);
      }
      a(i, j) = (a(i,j) - sum)/a(j,j);
    }
  }
}

/*!
  \brief
  Finds the next Newton iterate by the double-dogleg method.
   
  The new iterate is returned in the xpls vector. This method drives the
  dogstp function.

  Translated by Steve Verrill, April 15, 1998.

  \section template_args Template Parameters

  \param M  SCPPNT Matrix type.
  \param V  SCPPNT Vector type.
  \param FUNC Type of function to minimize.

  \section function_args Function Parameters

  \param[in]  &x         Address of old iterate (X[k-1]).
  \param[in]  &f         Address of function value at the old iterate.
  \param[in]  &g         Address of gradient or approximation at the old iterate.
  \param[in]  &a         Address of cholesky decomposition of hessian
      in lower triangular part and diagonal.
  \param[in]  &p         Address of newton step.
  \param[out] &xpls      Address of new iterate (X[k]).
  \param[out] &fpls      Address of function value at the new iterate.
  \param[in]  &sx        Address of scaling vector for x.
  \param[in,out] &dlt    Address of trust region radius (value needs to be retained
      between successive calls).
  \param[out] &iretcd    Address of return code: 0 -- satisfactory xpls found;
       1 -- failed to find satisfactory xpls sufficently distinct from x.
  \param[out] &mxtake    Address of boolean flag indicating that a step of maximum
      length was used length.
  \param[in]  &sc        Address of workspace (current step).
  \param[in]  &wrk1      Address of workspace (and place holding argument to tregup).
  \param[in]  &wrk2      Address of workspace.
  \param[in]  &wrk3      Address of workspace.
 */
template<class V, class M, class FUNC> void Uncmin<V, M, FUNC>::dogdrv(V &x, double &f, V &g, M &a,
    V &p, V &xpls, double &fpls, V &sx, double &dlt, int &iretcd, bool &mxtake, V &sc, V &wrk1,
    V &wrk2, V &wrk3)
{
  int i;
  double tmp, rnwtln;
  double fplsp;
  double cln;
  double eta;
  bool fstdog;
  bool nwtake;

  iretcd = 4;
  fstdog = true;
  tmp = 0.0;

  for (i = 1; i <= n; i++)
  {
    tmp += sx(i)*sx(i)*p(i)*p(i);
  }
  rnwtln = std::sqrt(tmp);

  while (iretcd > 1)
  {
    // FIND NEW STEP BY DOUBLE DOGLEG ALGORITHM
    dogstp(g, a, p, sx, rnwtln, dlt, nwtake, fstdog, wrk1, wrk2, cln, eta, sc);

    // CHECK NEW POINT AND UPDATE TRUST REGION
    tregup(x, f, g, a, sc, sx, nwtake, dlt, iretcd, wrk3, fplsp, xpls, fpls, mxtake, 2, wrk1);
  }
}

/*!
  \brief
  Finds the new step by the double-dogleg appproach.

  Translated by Steve Verrill, April 21, 1998.

  \section template_args Template Parameters

  \param M  SCPPNT Matrix type.
  \param V  SCPPNT Vector type.
  \param FUNC Type of function to minimize.

  \section function_args Function Parameters

  \param[in]      &g      Address of gradient at current iterate, g(X).
  \param[in]      &a      Address of cholesky decomposition of hessian in lower part.
      and diagonal.
  \param[in]      &p      Address of newton step.
  \param[in]      &sx     Address of scaling vector for x.
  \param[in]      rnwtln  Newton step length.
  \param[in,out] &dlt     Address of trust region radius.
  \param[in,out] &nwtake  Address of boolean, nwtake = TRUE if Newton step taken.
  \param[in,out] &fstdog  Address of boolean, fstdog = TRUE, if on first leg of dogleg.
  \param[in,out] &ssd     Address of workspace (Cauchy step to the minimum of the quadratic
      model in the scaled steeped descent direction; retaining value between 
      successive calls].
  \param[in,out] &v       Address of workspace (retainig value between successive calls).
  \param[in,out] &cln     Address of cauchy length (retaining value between successive calls).
  \param[in,out] &eta     Address of not documented, but "retains value between successive calls."
  \param[out]    &sc      Address of current step.
 */
template<class V, class M, class FUNC> void Uncmin<V, M, FUNC>::dogstp(V &g, M &a, V &p, V &sx,
    double rnwtln, double &dlt, bool &nwtake, bool &fstdog, V &ssd, V &v, double &cln, double &eta,
    V &sc)
{
  double alpha, beta, tmp, dot1, dot2, alam;
  int i, j;

  // CAN WE TAKE NEWTON STEP
  if (rnwtln <= dlt)
  {
    nwtake = true;

    for (i = 1; i <= n; i++)
    {
      sc(i) = p(i);
    }
    dlt = rnwtln;
  }
  else
  {
    // NEWTON STEP TOO LONG
    // CAUCHY STEP IS ON DOUBLE DOGLEG CURVE
    nwtake = false;

    if (fstdog)
    {
      // CALCULATE DOUBLE DOGLEG CURVE (SSD)
      fstdog = false;
      alpha = 0.0;

      for (i = 1; i <= n; i++)
      {
        alpha += (g(i)*g(i))/(sx(i)*sx(i));
      }

      beta = 0.0;

      for (i = 1; i <= n; i++)
      {
        tmp = 0.0;

        for (j = i; j <= n; j++)
        {
          tmp += (a(j, i)*g(j))/(sx(j)*sx(j));
        }
        beta += tmp*tmp;
      }

      for (i = 1; i <= n; i++)
      {
        ssd(i) = -(alpha/beta)*g(i)/sx(i);
      }

      cln = alpha*std::sqrt(alpha)/beta;

      eta = .2 + (.8*alpha*alpha)/(-beta*ddot(g, p));

      for (i = 1; i <= n; i++)
      {
        v(i) = eta*sx(i)*p(i) - ssd(i);
      }

      if (dlt == -1.0)
        dlt = min(cln, stepmx);
    }
    if (eta*rnwtln <= dlt)
    {
      // TAKE PARTIAL STEP IN NEWTON DIRECTION
      for (i = 1; i <= n; i++)
      {
        sc(i) = (dlt/rnwtln)*p(i);
      }
    }
    else
    {
      if (cln >= dlt)
      {
        // TAKE STEP IN STEEPEST DESCENT DIRECTION
        for (i = 1; i <= n; i++)
        {
          sc(i) = (dlt/cln)*ssd(i)/sx(i);
        }
      }
      else
      {
        // CALCULATE CONVEX COMBINATION OF SSD AND ETA*P
        // WHICH HAS SCALED LENGTH DLT
        dot1 = ddot(v, ssd); // Blas_f77.ddot_f77(n,v,1,ssd,1);
        dot2 = ddot(v, v); // Blas_f77.ddot_f77(n,v,1,v,1);

        alam = (-dot1 + std::sqrt((dot1*dot1) - dot2*(cln*cln - dlt*dlt)))/dot2;

        for (i = 1; i <= n; i++)
        {
          sc(i) = (ssd(i) + alam*v(i))/sx(i);
        }
      }
    }
  }
}

/*!
  \brief
  Solves Ax = b where A is a lower triangular matrix.

  If  is no longebr required by the calling routine, then the
  vectors b and x may share the same storage.

  Translated by Steve Verrill, April 21, 1998.

  \section template_args Template Parameters

  \param M  SCPPNT Matrix type.
  \param V  SCPPNT Vector type.
  \param FUNC Type of function to minimize.

  \section function_args Function Parameters

  \param[in]     &a     Address of the lower triangular matrix (preserved).
  \param[out]    &x     Address of solution vector.
  \param[in]     &b     Address of the right-hand side vector.
 */
template<class V, class M, class FUNC> void Uncmin<V, M, FUNC>::forslv(M &a, V &x, V &b)
{
  int i, im1, j;
  double sum;

  // SOLVE LX=B. (FORWARD SOLVE)
  x(1) = b(1)/a(1, 1);

  for (i = 2; i <= n; i++)
  {
    sum = 0.0;
    im1 = i - 1;

    for (j = 1; j <= im1; j++)
    {
      sum += a(i, j)*x(j);
    }
    x(i) = (b(i) - sum)/a(i,i);
  }
}

/*!
  \brief
  Finds a central difference approximation to the gradient g(x) of the function at point x.
 
  Translated by Steve Verrill, April 21, 1998.

  \section template_args Template Parameters

  \param M  SCPPNT Matrix type.
  \param V  SCPPNT Vector type.
  \param FUNC Type of function to minimize.

  \section function_args Function Parameters

  \param[in]     &x          Address of the point vector at which the gradient is to be approximated.
  \param[in]     &sx         Address of the scaling vector for x.
  \param[in]     rnoise      Relative noise in the function to be minimized.
  \param[out]    &g          Address of the central difference approximation to the gradient.
 */
template<class V, class M, class FUNC> void Uncmin<V, M, FUNC>::fstocd(V &x, V &sx, double rnoise,
    V &g)
{
  double stepi, xtempi, fplus, fminus, xmult;
  int i;

  // FIND I-TH  STEPSIZE, EVALUATE TWO NEIGHBORS IN DIRECTION OF I-TH
  // UNIT VECTOR, AND EVALUATE I-TH  COMPONENT OF GRADIENT.
  xmult = std::pow(rnoise, 1.0/3.0);

  for (i = 1; i <= n; i++)
  {
    stepi = xmult*max(std::fabs(x(i)), 1.0/sx(i));
    xtempi = x(i);

    x(i) = xtempi + stepi;
    fplus = minclass->f_to_minimize(x);

    x(i) = xtempi - stepi;
    fminus = minclass->f_to_minimize(x);

    x(i) = xtempi;

    g(i) = (fplus - fminus)/(2.0*stepi);
  }
}

/*!
  \brief
  Finds a next Newton iterate (xpls) by the More-Hebdon technique.
  
  This function drives hookst.
 
  Translated by Steve Verrill, April 23, 1998.

  \section template_args Template Parameters

  \param M  SCPPNT Matrix type.
  \param V  SCPPNT Vector type.
  \param FUNC Type of function to minimize.

  \section function_args Function Parameters

  \param[in]      &x          Address of old iterate (X[k-1]).
  \param[in]      &f          Address of function value at the old iterate.
  \param[in]      &g          Address of gradient or approximation at old iterate.
  \param[in]      &a          Address of cholesky decomposition of hessian in lower 
      triangle and diagonal. Upper triangle contains the off-diagonal elements of 
      the hessian (diagonal of hessian expected in udiag).
  \param[in]      &udiag      Address of vector containing diagonal of hessian.
  \param[in]      &p          Address of newton step.
  \param[out]     &xpls       Address of new iterate (X[k]).
  \param[out]     &fpls       Address of function value at the new iterate.
  \param[in]      &sx         Address of scaling vector for x.
  \param[in,out]  &dlt        Address of trust region radius [data member].
  \param[out]     &iretcd     Address of return code: iretcd = 0, satisfactory xpls found;
      iretcd = 1, failed to find satisfactory xpls sufficiently distinct from x.
  \param[out]     &mxtake     Address of boolean flag indicating step of maximum length used.
  \param[in,out]  &amu        Address of ? [Retaining value between successive calls].
  \param[in,out]  &dltp       Address of ? [Retaining value between successive calls].
  \param[in,out]  &phi        Address of ? [Retaining value between successive calls].
  \param[in,out]  &phip0      Address of ? [Retaining value between successive calls].
  \param[in]      &sc         Address of workspace.
  \param[in]      &xplsp      Address of workspace.
  \param[in]      &wrk0       Address of workspace.
  \param[in]      &itncnt     Address of iteration count.
 */
template<class V, class M, class FUNC> void Uncmin<V, M, FUNC>::hookdr(V &x, double &f, V &g, M &a,
    V &udiag, V &p, V &xpls, double &fpls, V &sx, double &dlt, int &iretcd, bool &mxtake,
    double &amu, double &dltp, double &phi, double &phip0, V &sc, V &xplsp, V &wrk0, int &itncnt)
{
  int i, j;
  bool fstime;
  bool nwtake;
  double tmp, rnwtln, alpha, beta;

  double fplsp;

  iretcd = 4;
  fstime = true;

  tmp = 0.0;

  for (i = 1; i <= n; i++)
  {
    tmp += sx(i)*sx(i)*p(i)*p(i);
  }
  rnwtln = std::sqrt(tmp);

  if (itncnt == 1)
  {
    amu = 0.0;

    // IF FIRST ITERATION AND TRUST REGION NOT PROVIDED BY USER,
    // COMPUTE INITIAL TRUST REGION.
    if (dlt == -1.0)
    {
      alpha = 0.0;

      for (i = 1; i <= n; i++)
      {
        alpha += (g(i)*g(i))/(sx(i)*sx(i));
      }
      beta = 0.0;

      for (i = 1; i <= n; i++)
      {
        tmp = 0.0;

        for (j = i; j <= n; j++)
        {
          tmp += (a(j, i)*g(j))/(sx(j)*sx(j));
        }
        beta += tmp*tmp;
      }
      dlt = alpha*std::sqrt(alpha)/beta;
      dlt = min(dlt, stepmx);
    }
  }
  while (iretcd > 1)
  {
    // FIND NEW STEP BY MORE-HEBDON ALGORITHM
    hookst(g, a, udiag, p, sx, rnwtln, dlt, amu, dltp, phi, phip0, fstime, sc, nwtake, wrk0);

    dltp = dlt;

    // CHECK NEW POINT AND UPDATE TRUST REGION
    tregup(x, f, g, a, sc, sx, nwtake, dlt, iretcd, xplsp, fplsp, xpls, fpls, mxtake, 3, udiag);
  }
}

/*!
  \brief
  Finds a new step (xpls) by the More-Hebdon algorithm.
  
  This function is driven by hookdr.
  
  Translated by Steve Verrill, April 24, 1998.

  \section template_args Template Parameters

  \param M  SCPPNT Matrix type.
  \param V  SCPPNT Vector type.
  \param FUNC Type of function to minimize.

  \section function_args Function Parameters

  \param[in]      &g          Address of gradient at the current iterate, g(X).
  \param[in]      &a          Address of cholesky decomposition of hessian in lower 
      triangle and diagonal. Upper triangle contains the off-diagonal elements of 
      the hessian (diagonal of hessian expected in udiag).
  \param[in]      &udiag      Address of vector containing diagonal of hessian.
  \param[in]      &p          Address of newton step.
  \param[in]      &sx         Address of scaling vector for x.
  \param[in]      rnwtln      Newton step length.
  \param[in,out]  &dlt        Address of trust region radius at last exit from this routine.
  \param[in,out]  &amu        Address of ? [Retaining value between successive calls].
  \param[in]      &dltp       Address of trust region radius at last exit from this routine.
  \param[in,out]  &phi        Address of ? [Retaining value between successive calls].
  \param[in,out]  &phip0      Address of ? [Retaining value between successive calls].
  \param[in,out]  &fstime     Address of boolean flag: fstime = TRUE, if first entry to this 
      routine during the k-th iteration.
  \param[out]     &sc         Address of current step.
  \param[out]     &nwtake     Address of boolean flag: nwtake = TRUE, if newton step taken.
  \param[in]      &wrk0       Address of workspace.
 */
template<class V, class M, class FUNC> void Uncmin<V, M, FUNC>::hookst(V &g, M &a, V &udiag, V &p,
    V &sx, double rnwtln, double &dlt, double &amu, double &dltp, double &phi, double &phip0,
    bool &fstime, V &sc, bool &nwtake, V &wrk0)
{
  double hi, alo;
  double phip, amulo, amuup, stepln;
  int i, j;
  bool done;

  double addmax;

  // HI AND ALO ARE CONSTANTS USED IN THIS ROUTINE.
  // CHANGE HERE IF OTHER VALUES ARE TO BE SUBSTITUTED.
  phip = 0.0;
  hi = 1.5;
  alo = .75;

  if (rnwtln <= hi*dlt)
  {
    //       TAKE NEWTON STEP
    nwtake = true;

    for (i = 1; i <= n; i++)
    {
      sc(i) = p(i);
    }

    dlt = min(dlt, rnwtln);
    amu = 0.0;

    return;
  }
  else
  {
    // NEWTON STEP NOT TAKEN
    nwtake = false;

    if (amu > 0.0)
    {
      amu -= (phi + dltp)*((dltp - dlt) + phi)/(dlt*phip);
    }

    phi = rnwtln - dlt;

    if (fstime)
    {
      for (i = 1; i <= n; i++)
      {
        wrk0(i) = sx(i)*sx(i)*p(i);
      }
      // SOLVE L*Y = (SX**2)*P
      forslv(a, wrk0, wrk0);
      phip0 = -std::pow(dnrm2(wrk0), 2)/rnwtln;
      fstime = false;
    }
    phip = phip0;
    amulo = -phi/phip;

    amuup = 0.0;

    for (i = 1; i <= n; i++)
    {
      amuup += (g(i)*g(i))/(sx(i)*sx(i));
    }
    amuup = std::sqrt(amuup)/dlt;

    done = false;

    // TEST VALUE OF AMU; GENERATE NEXT AMU IF NECESSARY
    while (!done)
    {
      if (amu < amulo || amu > amuup)
      {
        amu = max(std::sqrt(amulo*amuup),amuup*.001);
      }

      // COPY (H,UDIAG) TO L
      // WHERE H <-- H+AMU*(SX**2) [DO NOT ACTUALLY CHANGE (H,UDIAG)]
      for (j = 1; j <= n; j++)
      {
        a(j,j) = udiag(j) + amu*sx(j)*sx(j);

        for (i = j + 1; i <= n; i++)
        {
          a(i,j) = a(j,i);
        }
      }

      // FACTOR H=L(L+)
      choldc(a,0.0,std::sqrt(epsm),addmax);

      // SOLVE H*P = L(L+)*SC = -G
      for (i = 1; i <= n; i++)
      {
        wrk0(i) = -g(i);
      }

      lltslv(a,sc,wrk0);

      // RESET H.  NOTE SINCE UDIAG HAS NOT BEEN DESTROYED WE NEED DO
      // NOTHING HERE.  H IS IN THE UPPER PART AND IN UDIAG, STILL INTACT
      stepln = 0.0;

      for (i = 1; i <= n; i++)
      {
        stepln += sx(i)*sx(i)*sc(i)*sc(i);
      }

      stepln = std::sqrt(stepln);
      phi = stepln - dlt;

      for (i = 1; i <= n; i++)
      {
        wrk0(i) = sx(i)*sx(i)*sc(i);
      }

      forslv(a,wrk0,wrk0);

      phip = -std::pow(dnrm2(wrk0),2)/stepln;

      if ((alo*dlt <= stepln && stepln <= hi*dlt) ||
      (amuup-amulo <= 0.0))
      {
        // SC IS ACCEPTABLE HOOKSTEP
        done = true;
      }
      else
      {
        // SC NOT ACCEPTABLE HOOKSTEP.  SELECT NEW AMU
        amulo = max(amulo,amu-(phi/phip));
        if (phi < 0.0) amuup = min(amuup,amu);
        amu -= (stepln*phi)/(dlt*phip);
      }
    }
  }
}

/*!
  \brief
  Solves Ax = b where A has the form LL', and only the lower triangular part, L, is stored.

  Note: If b is not required by the calling program, then b and x may share the same storage.
 
  Translated by Steve Verrill, April 27, 1998.

  \section template_args Template Parameters

  \param M  SCPPNT Matrix type.
  \param V  SCPPNT Vector type.
  \param FUNC Type of function to minimize.

  \section function_args Function Parameters

  \param[in]      &a     Address of matrix of form L (A= LL'). On return, this matrix is unchanged.
  \param[out]     &x     Address of solution vector.
  \param[in]      &b     Address of the right-hand side vector.
 */
template<class V, class M, class FUNC> void Uncmin<V, M, FUNC>::lltslv(M &a, V &x, V &b)
{
  // FORWARD SOLVE, RESULT IN X
  forslv(a, x, b);
  // BACK SOLVE, RESULT IN X
  bakslv(a, x, x);
}

/*!
  \brief
  Finds a next Newton iterate by line search.
  
  Translated by Steve Verrill, May 15, 1998.

  \section template_args Template Parameters

  \param M  SCPPNT Matrix type.
  \param V  SCPPNT Vector type.
  \param FUNC Type of function to minimize.

  \section function_args Function Parameters

  \param[in]  &x      Address of old iterate, (X[k-1]).
  \param[in]  &f      Address of function value at old iterate.
  \param[in]  &g      Address of gradient or approximation at old iterate.
  \param[in]  &p      Address of non-zero Newton step.
  \param[out] &xpls   Address of new iterate. (X[k]).
  \param[out] &fpls   Address of function value at new iterate.
  \param[out] &mxtake Address of boolean flag indicating whether the step of
      maximum length was used.
  \param[out] &iretcd Address of return code
  \param[in]  &sx     Address of scaling vector for x.
 */
template<class V, class M, class FUNC> void Uncmin<V, M, FUNC>::lnsrch(V &x, double &f, V &g, V &p,
    V &xpls, double &fpls, bool &mxtake, int &iretcd, V &sx)
{
  int i;
  double tmp, sln, scl, slp, rln, rmnlmb, almbda, tlmbda;
  double t1, t2, t3, a, b, disc, pfpls, plmbda;

  pfpls = 0.0;
  plmbda = 0.0;

  mxtake = false;
  iretcd = 2;

  tmp = 0.0;

  for (i = 1; i <= n; i++)
  {
    tmp += sx(i)*sx(i)*p(i)*p(i);
  }
  sln = std::sqrt(tmp);

  if (sln > stepmx)
  {
    // NEWTON STEP LONGER THAN MAXIMUM ALLOWED
    scl = stepmx/sln;
    p *= scl; // Uncmin_f77.sclmul_f77(n,scl,p,p);
    sln = stepmx;
  }
  slp = ddot(g, p);
  rln = 0.0;

  for (i = 1; i <= n; i++)
  {
    rln = max(rln, std::fabs(p(i))/max(std::fabs(x(i)), 1.0/sx(i)));
  }

  rmnlmb = steptl/rln;
  almbda = 1.0;

  // LOOP
  // CHECK IF NEW ITERATE SATISFACTORY.  GENERATE NEW LAMBDA IF NECESSARY.
  int iteration = 0; // BAH
  while (iretcd >= 2)
  {
    for (i = 1; i <= n; i++)
    {
      xpls(i) = x(i) + almbda*p(i);
    }
    /* Find valid function parameters - BAH */
    if (!(minclass->ValidParameters(xpls)))
    {
      almbda *= 0.1;
      if (almbda < rmnlmb)
      {
        iretcd = 1; // no satisfactory xpls found sufficiently distinct from x
        return;
      }
      continue;
    }

    fpls = minclass->f_to_minimize(xpls);

    if (fpls <= (f + slp*.0001*almbda))
    {
      // SOLUTION FOUND
      iretcd = 0;
      if (almbda == 1.0 && sln > .99*stepmx)
        mxtake = true;
    }
    else
    {
      // SOLUTION NOT (YET) FOUND
      // test for iteration limit to prevent infinate loop when
      // fpls of new iterate is NaN - BAH
      if (++iteration > 100)
      {
        iretcd = -7;
        if (fPrintResults)
        {
          std::fprintf(mfile, "\n\n\nLNSRCH      Number of iterations exceeded.\n");
          std::fprintf(mfile, "LNSRCH      fpls = %e\n", fpls);
        }
        return;
      } // BAH

      if (almbda < rmnlmb)
      {
        // NO SATISFACTORY XPLS FOUND SUFFICIENTLY DISTINCT FROM X
        iretcd = 1;
      }
      else
      {
        // CALCULATE NEW LAMBDA
        if (iteration == 1)
        {
          // FIRST BACKTRACK: QUADRATIC FIT
          tlmbda = -slp/(2.0*(fpls - f - slp));
        }
        else
        {
          // ALL SUBSEQUENT BACKTRACKS: CUBIC FIT
          t1 = fpls - f - almbda*slp;
          t2 = pfpls - f - plmbda*slp;
          t3 = 1.0/(almbda - plmbda);
          a = t3*(t1/(almbda*almbda) - t2/(plmbda*plmbda));
          b = t3*(t2*almbda/(plmbda*plmbda) - t1*plmbda/(almbda*almbda));
          disc = b*b - 3.0*a*slp;

          if (disc > b*b)
          {
            // ONLY ONE POSITIVE CRITICAL POINT, MUST BE MINIMUM
            double signone = (a < 0.0) ? -1.0 : 1.0;
            tlmbda = (-b + signone*std::sqrt(disc))/(3.0*a);
          }
          else
          {
            // BOTH CRITICAL POINTS POSITIVE, FIRST IS MINIMUM
            double signone = (a < 0.0) ? -1.0 : 1.0;
            tlmbda = (-b - signone*std::sqrt(disc))/(3.0*a);
          }

          if (tlmbda > .5*almbda)
            tlmbda = .5*almbda;
        }
        plmbda = almbda;
        pfpls = fpls;

        if (tlmbda < almbda/10.0)
        {
          almbda *= .1;
        }
        else
        {
          almbda = tlmbda;
        }
      }
    }
  }
}

// newton direction
template<class V, class M, class FUNC> void Uncmin<V, M, FUNC>::lnsrch(V &x, double &f, V &g, V &p,
	V &xpls, double &fpls, bool &mxtake, int &iretcd, V &sx, int normal_x, int normal_y, int normal_z)
{
	int i;
	double tmp, sln, scl, slp, rln, rmnlmb, almbda, tlmbda;
	double t1, t2, t3, a, b, disc, pfpls, plmbda;

	pfpls = 0.0;
	plmbda = 0.0;

	mxtake = false;
	iretcd = 2;

	tmp = 0.0;

	for (i = 1; i <= n; i++)
	{
		tmp += sx(i)*sx(i)*p(i)*p(i);
	}
	sln = std::sqrt(tmp);

	if (sln > stepmx)
	{
		// NEWTON STEP LONGER THAN MAXIMUM ALLOWED
		scl = stepmx / sln;
		p *= scl; // Uncmin_f77.sclmul_f77(n,scl,p,p);
		sln = stepmx;
	}
	slp = ddot(g, p);//g only used here
	rln = 0.0;

	for (i = 1; i <= n; i++)
	{
		rln = max(rln, std::fabs(p(i)) / max(std::fabs(x(i)), 1.0 / sx(i)));
	}

	rmnlmb = steptl / rln;
	almbda = 1.0;

	// LOOP
	// CHECK IF NEW ITERATE SATISFACTORY.  GENERATE NEW LAMBDA IF NECESSARY.
	int iteration = 0; // BAH
	while (iretcd >= 2)
	{
		for (i = 1; i <= n; i++)
		{
			xpls(i) = x(i) + almbda * p(i);
			//std::cout << "xpls:" << i << " : " << x(i) << std::endl;
		}
		/* Find valid function parameters - BAH */
		if (!(minclass->ValidParameters(xpls)))
		{
			almbda *= 0.1;
			if (almbda < rmnlmb)
			{
				iretcd = 1; // no satisfactory xpls found sufficiently distinct from x
				return;
			}
			continue;
		}

		//fpls = minclass->f_to_minimize(xpls);
		fpls = minclass->f_to_minimize(xpls, normal_x, normal_y, normal_z);

		if (fpls <= (f + slp * .0001*almbda))
		{
			// SOLUTION FOUND
			iretcd = 0;
			if (almbda == 1.0 && sln > .99*stepmx)
				mxtake = true;
		}
		else
		{
			// SOLUTION NOT (YET) FOUND
			// test for iteration limit to prevent infinate loop when
			// fpls of new iterate is NaN - BAH
			if (++iteration > 100)
			{
				iretcd = -7;
				if (fPrintResults)
				{
					std::fprintf(mfile, "\n\n\nLNSRCH      Number of iterations exceeded.\n");
					std::fprintf(mfile, "LNSRCH      fpls = %e\n", fpls);
				}
				return;
			} // BAH

			if (almbda < rmnlmb)
			{
				// NO SATISFACTORY XPLS FOUND SUFFICIENTLY DISTINCT FROM X
				iretcd = 1;
			}
			else
			{
				// CALCULATE NEW LAMBDA
				if (iteration == 1)
				{
					// FIRST BACKTRACK: QUADRATIC FIT
					tlmbda = -slp / (2.0*(fpls - f - slp));
				}
				else
				{
					// ALL SUBSEQUENT BACKTRACKS: CUBIC FIT
					t1 = fpls - f - almbda * slp;
					t2 = pfpls - f - plmbda * slp;
					t3 = 1.0 / (almbda - plmbda);
					a = t3 * (t1 / (almbda*almbda) - t2 / (plmbda*plmbda));
					b = t3 * (t2*almbda / (plmbda*plmbda) - t1 * plmbda / (almbda*almbda));
					disc = b * b - 3.0*a*slp;

					if (disc > b*b)
					{
						// ONLY ONE POSITIVE CRITICAL POINT, MUST BE MINIMUM
						double signone = (a < 0.0) ? -1.0 : 1.0;
						
						if (disc < 0)
						{
							iretcd = -7;
							return;
						}
						
						tlmbda = (-b + signone * std::sqrt(disc)) / (3.0*a);
					}
					else
					{
						// BOTH CRITICAL POINTS POSITIVE, FIRST IS MINIMUM
						double signone = (a < 0.0) ? -1.0 : 1.0;

						if (disc < 0)
						{
							iretcd = -7;
							return;
						}
						tlmbda = (-b - signone * std::sqrt(disc)) / (3.0*a);
					}

					if (tlmbda > .5*almbda)
						tlmbda = .5*almbda;
				}
				plmbda = almbda;
				pfpls = fpls;

				if (tlmbda < almbda / 10.0)
				{
					almbda *= .1;
				}
				else
				{
					almbda = tlmbda;
				}
			}
		}
	}
}
/*!
  \brief
  Computes y = Lx where L is a lower triangular matrix stored in A.
 
  Translated by Steve Verrill, April 27, 1998.

  \section template_args Template Parameters

  \param M  SCPPNT Matrix type.
  \param V  SCPPNT Vector type.
  \param FUNC Type of function to minimize.

  \section function_args Function Parameters

  \param[in]      &a     Address of lower triangular matrix A (containing L).
  \param[in]      &x     Address of operand vector.
  \param[out]     &y     Address of result vector.
 */
template<class V, class M, class FUNC> void Uncmin<V, M, FUNC>::mvmltl(M &a, V &x, V &y)
{
  double sum;
  int i, j;

  for (i = 1; i <= n; i++)
  {
    sum = 0.0;
    for (j = 1; j <= i; j++)
    {
      sum += a(i, j)*x(j);
    }
    y(i) = sum;
  }
}

/*!
  \brief
  Computes y = Ax where A is a symmetric matrix stored in its lower triangular part; x and y are vectors.
 
  Translated by Steve Verrill, April 27, 1998.

  \section template_args Template Parameters

  \param M  SCPPNT Matrix type.
  \param V  SCPPNT Vector type.
  \param FUNC Type of function to minimize.

  \section function_args Function Parameters

  \param[in]      &a     Address of the symmetric matrix A.
  \param[in]      &x     Address of operand vector x.
  \param[out]     &y     Address of result vector y.
 */
template<class V, class M, class FUNC> void Uncmin<V, M, FUNC>::mvmlts(M &a, V &x, V &y)
{
  double sum;
  int i, j;

  for (i = 1; i <= n; i++)
  {
    sum = 0.0;

    for (j = 1; j <= i; j++)
    {
      sum += a(i, j)*x(j);
    }
    for (j = i + 1; j <= n; j++)
    {
      sum += a(j, i)*x(j);
    }
    y(i) = sum;
  }
}

/*!
  \brief
  Computes y = L'x where L is a lower triangular matrix stored in A; y and x are vectors.
 
  Note:  The L transpose is taken implicitly.
  
  Translated by Steve Verrill, April 27, 1998.

  \section template_args Template Parameters

  \param M  SCPPNT Matrix type.
  \param V  SCPPNT Vector type.
  \param FUNC Type of function to minimize.

  \section function_args Function Parameters

  \param[in]      &a     Address of the lower triangular matrix A.
  \param[in]      &x     Address of operand vector x.
  \param[out]     &y     Address of result vector y.
 */
template<class V, class M, class FUNC> void Uncmin<V, M, FUNC>::mvmltu(M &a, V &x, V &y)
{
  double sum;
  int i, j;

  for (i = 1; i <= n; i++)
  {
    sum = 0.0;
    for (j = i; j <= n; j++)
    {
      sum += a(j, i)*x(j);
    }
    y(i) = sum;
  }
}

/*!
  \brief
  Interchanges rows i and i+1 of the upper Hessenberg matrix r, in columns i to n.

  Translated by Steve Verrill, April 29, 1998.

  \section template_args Template Parameters

  \param M  SCPPNT Matrix type.
  \param V  SCPPNT Vector type.
  \param FUNC Type of function to minimize.

  \section function_args Function Parameters

  \param[in,out]  &r     Address of upper Hessenberg matrix r.
  \param[in]      i      Index of row to interchange (i < n)
 */
template<class V, class M, class FUNC> void Uncmin<V, M, FUNC>::qraux1(M &r, int i)
{
  double tmp;
  int j, ip1;

  ip1 = i + 1;

  for (j = i; j <= n; j++)
  {
    tmp = r(i, j);
    r(i, j) = r(ip1,j);
    r(ip1, j) = tmp;
  }
}

/*!
  \brief
  Pre-multiplies the upper Hessenberg matrix r by the Jacobi rotation j(i,i+1,a,b).
  
  Translated by Steve Verrill, April 29, 1998.

  \section template_args Template Parameters

  \param M  SCPPNT Matrix type.
  \param V  SCPPNT Vector type.
  \param FUNC Type of function to minimize.

  \section function_args Function Parameters

  \param[in,out]  &r    Address of upper Hessenberg matrix r.
  \param[in]      i     Index of row.
  \param[in]      a     scalar.
  \param[in]      b     scalar.
 */
template<class V, class M, class FUNC> void Uncmin<V, M, FUNC>::qraux2(M &r, int i, double a,
    double b)
{
  double den, c, s, y, z;
  int j, ip1;

  ip1 = i + 1;

  den = std::sqrt(a*a + b*b);
  c = a/den;
  s = b/den;

  for (j = i; j <= n; j++)
  {
    y = r(i, j);
    z = r(ip1, j);
    r(i, j) = c*y - s*z;
    r(ip1, j) = s*y + c*z;
  }
}


/*!
  \brief
  Finds an orthogonal n by n matrix, Q*, and an upper triangular n by n matrix, R*, such that
  (Q*)(R*) = R+u*v'.
  
  Translated by Steve Verrill, May 11, 1998.

  \section template_args Template Parameters

  \param M  SCPPNT Matrix type.
  \param V  SCPPNT Vector type.
  \param FUNC Type of function to minimize.

  \section function_args Function Parameters

  \param[in,out]  &a    Address of matrix A. Contains R on input; contains R* at output.
  \param[in]      &u    Address of input vector u.
  \param[in]      &v    Address of input vector v.
 */
template<class V, class M, class FUNC> void Uncmin<V, M, FUNC>::qrupdt(M &a, V &u, V &v)
{
  int k, km1, ii, i, j;
  double t1, t2;

  // DETERMINE LAST NON-ZERO IN U(.)
  k = n;

  while (u(k) == 0 && k > 1)
  {
    k--;
  }

  // (K-1) JACOBI ROTATIONS TRANSFORM
  // R + U(V+) --> (R*) + (U(1)*E1)(V+)
  // WHICH IS UPPER HESSENBERG
  km1 = k - 1;

  for (ii = 1; ii <= km1; ii++)
  {
    i = km1 - ii + 1;

    if (u(i) == 0.0)
    {
      qraux1(a, i);
      u(i) = u(i+1);
    }
    else
    {
      qraux2(a, i, u(i), -u(i+1));
      u(i) = std::sqrt(u(i)*u(i) + u(i+1)*u(i+1));
    }
  }

  // R <-- R + (U(1)*E1)(V+)
  for (j = 1; j <= n; j++)
  {
    a(1, j) += u(1)*v(j);
  }

  // (K-1) JACOBI ROTATIONS TRANSFORM UPPER HESSENBERG R
  // TO UPPER TRIANGULAR (R*)
  km1 = k - 1;

  for (i = 1; i <= km1; i++)
  {
    if (a(i, i) == 0.0)
    {
      qraux1(a, i);
    }
    else
    {
      t1 = a(i, i);
      t2 = -a(i+1, i);

      qraux2(a, i, t1, t2);
    }
  }
}

/*!
  \brief
  Updates the Hessian by the BFGS factored technique.
  
  Translated by Steve Verrill, May 14, 1998.

  \section template_args Template Parameters

  \param M  SCPPNT Matrix type.
  \param V  SCPPNT Vector type.
  \param FUNC Type of function to minimize.

  \section function_args Function Parameters

  \param[in]      &x       Address of old iterate (X[k-1]).
  \param[in]      &g       Address of gradient or approximation at the old iterate.
  \param[in,out]  &a       Address of matrix A. A contains on entry the cholesky decomposition 
      of the hessian in its lower triangle and diagonal; at exit, A contains the updated 
      cholesky decomposition of the hessian in its lower triangle and diagonal.
  \param[in]      &xpls    Address of the new iterate (X[k]).
  \param[in]      &gpls    Address of gradient or approximation at the new iterate.
  \param[in]      &itncnt  Address of variable containing iteration count.
  \param[in]      rnf      Relative noise in optimization function f_to_minimize.
  \param[in]      &iagflg  Address of integer Flag: iagflag = 1 if an analytic gradient is 
      supplied, 0 otherwise.
  \param[in,out]  &noupdt  Address of boolean flag: noupdt = TRUE, when no update yet (retain 
      value between successive calls).
  \param[in]      &s       Address of Workspace vector s.
  \param[in]      &y       Address of Workspace vector y.
  \param[in]      &u       Address of Workspace vector u.
  \param[in]      &w       Address of Workspace vector w.
 */
template<class V, class M, class FUNC> void Uncmin<V, M, FUNC>::secfac(V &x, V &g, M &a, V &xpls,
    V &gpls, int &itncnt, double rnf, int &iagflg, bool &noupdt, V &s, V &y, V &u, V &w)
{
  bool skpupd;
  int i, j, im1;
  double den1, snorm2, ynrm2, den2, alp, reltol;

  if (itncnt == 1)
    noupdt = true;

  for (i = 1; i <= n; i++)
  {
    s(i) = xpls(i) - x(i);
    y(i) = gpls(i) - g(i);
  }

  den1 = ddot(s, y); // Blas_f77.ddot_f77(n,s,1,y,1);
  snorm2 = dnrm2(s); // Blas_f77.dnrm2_f77(n,s,1);
  ynrm2 = dnrm2(y); // Blas_f77.dnrm2_f77(n,y,1);

  if (den1 >= std::sqrt(epsm)*snorm2*ynrm2)
  {
    mvmltu(a, s, u);
    den2 = ddot(u, u); // Blas_f77.ddot_f77(n,u,1,u,1);

    // L <-- SQRT(DEN1/DEN2)*L
    alp = std::sqrt(den1/den2);

    if (noupdt)
    {
      for (j = 1; j <= n; j++)
      {
        u(j) *= alp;

        for (i = j; i <= n; i++)
        {
          a(i, j) *= alp;
        }
      }
      noupdt = false;
      // den2 = den1;  This value of den2 is never used - BAH
      alp = 1.0;
    }
    skpupd = true;

    // W = L(L+)S = HS
    mvmltl(a, u, w);

    i = 1;

    if (iagflg == 0)
    {
      reltol = std::sqrt(rnf);
    }
    else
    {
      reltol = rnf;
    }

    while (i <= n && skpupd)
    {
      if (std::fabs(y(i) - w(i)) >= reltol*max(std::fabs(g(i)), std::fabs(gpls(i))))
      {
        skpupd = false;
      }
      else
      {
        i++;
      }
    }
    if (!skpupd)
    {
      // W=Y-ALP*L(L+)S
      for (i = 1; i <= n; i++)
      {
        w(i) = y(i) - alp*w(i);
      }
      // ALP=1/SQRT(DEN1*DEN2)
      alp /= den1;

      // U=(L+)/SQRT(DEN1*DEN2) = (L+)S/SQRT((Y+)S * (S+)L(L+)S)
      for (i = 1; i <= n; i++)
      {
        u(i) *= alp;
      }

      // COPY L INTO UPPER TRIANGULAR PART.  ZERO L.
      for (i = 2; i <= n; i++)
      {
        im1 = i - 1;

        for (j = 1; j <= im1; j++)
        {
          a(j, i) = a(i,j);
          a(i, j) = 0.0;
        }
      }

      // FIND Q, (L+) SUCH THAT  Q(L+) = (L+) + U(W+)
      qrupdt(a, u, w);

      // UPPER TRIANGULAR PART AND DIAGONAL OF A NOW CONTAIN UPDATED
      // CHOLESKY DECOMPOSITION OF HESSIAN.  COPY BACK TO LOWER
      // TRIANGULAR PART.
      for (i = 2; i <= n; i++)
      {
        im1 = i - 1;

        for (j = 1; j <= im1; j++)
        {
          a(i, j) = a(j,i);
        }
      }
    }
  }
}


/*!
  \brief
  Updates the Hessian by the BFGS unfactored approach.
  
  Translated by Steve Verrill, May 8, 1998.

  \section template_args Template Parameters

  \param M  SCPPNT Matrix type.
  \param V  SCPPNT Vector type.
  \param FUNC Type of function to minimize.

  \section function_args Function Parameters

  \param[in]      &x       Address of old iterate (X[k-1]).
  \param[in]      &g       Address of gradient or approximation at the old iterate.
  \param[in,out]  &a       Address of matrix A. On entry, A contains the off-diagonal entries
      of the approximate Hhssian at the old iterate in its upper triangle (Note: The diagonal 
      of the hessian is passed in the udiag vector). At exit, A contains the updated approximate
      hessian at the new iterate in its lower triangle and diagonal.
  \param[in]      &udiag   Address of vector containing diagonal of Hessian at old iterate.    
  \param[in]      &xpls    Address of the new iterate (X[k]).
  \param[in]      &gpls    Address of gradient or approximation at the new iterate.
  \param[in]      &itncnt  Address of variable containing iteration count.
  \param[in]      rnf      Relative noise in optimization function f_to_minimize.
  \param[in]      &iagflg  Address of integer Flag: iagflag = 1 if an analytic gradient is 
      supplied, 0 otherwise.
  \param[in,out]  &noupdt  Address of boolean flag: noupdt = TRUE, when no update yet (retain 
      value between successive calls).
  \param[in]      &s       Address of Workspace vector s.
  \param[in]      &y       Address of Workspace vector y.
  \param[in]      &t       Address of Workspace vector t.
 */
template<class V, class M, class FUNC> void Uncmin<V, M, FUNC>::secunf(V &x, V &g, M &a, V &udiag,
    V &xpls, V &gpls, int &itncnt, double rnf, int &iagflg, bool &noupdt, V &s, V &y, V &t)
{
  double den1, snorm2, ynrm2, den2, gam, tol;
  int i, j;
  bool skpupd;

  // COPY HESSIAN IN UPPER TRIANGULAR PART AND UDIAG TO
  // LOWER TRIANGULAR PART AND DIAGONAL
  for (j = 1; j <= n; j++)
  {
    a(j, j) = udiag(j);

    for (i = j+1; i <= n; i++)
    {
      a(i, j) = a(j,i);
    }
  }
  if (itncnt == 1)
    noupdt = true;

  for (i = 1; i <= n; i++)
  {
    s(i) = xpls(i) - x(i);
    y(i) = gpls(i) - g(i);
  }

  den1 = ddot(s, y); // Blas_f77.ddot_f77(n,s,1,y,1);

  snorm2 = dnrm2(s); // Blas_f77.dnrm2_f77(n,s,1);

  ynrm2 = dnrm2(y); // Blas_f77.dnrm2_f77(n,y,1);

  if (den1 >= std::sqrt(epsm)*snorm2*ynrm2)
  {
    mvmlts(a, s, t);

    den2 = ddot(s, t); // Blas_f77.ddot_f77(n,s,1,t,1);

    if (noupdt)
    {
      // H <-- [(S+)Y/(S+)HS]H
      gam = den1/den2;

      den2 = gam*den2;

      for (j = 1; j <= n; j++)
      {
        t(j) *= gam;

        for (i = j; i <= n; i++)
        {
          a(i, j) *= gam;
        }
      }
      noupdt = false;
    }

    skpupd = true;

    // CHECK UPDATE CONDITION ON ROW I
    for (i = 1; i <= n; i++)
    {
      tol = rnf*max(std::fabs(g(i)), std::fabs(gpls(i)));
      if (iagflg == 0)
        tol /= std::sqrt(rnf);

      if (std::fabs(y(i) - t(i)) >= tol)
      {
        skpupd = false;
        break;
      }
    }

    if (!skpupd)
    {
      // BFGS UPDATE
      for (j = 1; j <= n; j++)
      {
        for (i = j; i <= n; i++)
        {
          a(i, j) += y(i)*y(j)/den1 - t(i)*t(j)/den2;
        }
      }
    }
  }
}

/*!
  \brief
  Determines whether to accept xpls = x + sc as the next iterate and update the trust region dlt.

  Translated by Steve Verrill, May 11, 1998.

  \section template_args Template Parameters

  \param M  SCPPNT Matrix type.
  \param V  SCPPNT Vector type.
  \param FUNC Type of function to minimize.

  \section function_args Function Parameters

  \param[in]      &x        Address of old iterate (X[k-1]).
  \param[in]      &f        Address of Function value at old iterate, f(X).
  \param[in]      &g        Address of Gradient or approximation at old iterate, G(X).
  \param[in]      &a        Address of matrix A, containing cholesky decomposition of 
      hessian in lower triangular part and diagonal. The upper triangular section of A
      contains the off-diagonal entries of the hessian itself.
  \param[in]      &sc       Address of the current step.
  \param[in]      &sx       Address of the scaling vector for x.
  \param[in]      &nwtake   Address of boolean flag: nwtake = TRUE, if Newton step taken.
  \param[in,out]  &dlt      Address of trust region radius. 
  \param[in,out]  &iretcd   Address of return code. Explanation of iretcd values are:
      If iretcd = 0, then xpls is accepted as next iterate, dlt is the trust region 
      radius for the next iteration;
      If iretcd = 1, then xpls is unsatisfactory but accepted as next iterate because 
      xpls - x is less than the smallest allowable step length;
      If iretcd = 2, then f(xpls) is too large -- the current iteration is continued with new 
      reduced dlt;
      If iretcd = 3, then f(xpls) sufficiently small, but the quadratic model predicts f(xpls) 
      sufficiently well to continue current iteration with a new doubled dlt.
  \param[in,out]  &xplsp    Address of Workspace (value needs to be retained between
 *                  successive calls of k-th global step).
  \param[in,out]  &fplsp    Address of ? [Retaining between successive calls].
  \param[out]     &xpls     Address of vector containing new iterate (X[k]).
  \param[out]     &fpls     Address of Function value at new iterate.
  \param[out]     &mxtake   Address of Boolean flag indicating step of maximum length used
  \param[in]      method    Algorithm to use to solve minimization problem. Values are:
      method = 1, line search;
      method = 2, double dogleg;
      methos = 3, More-Hebdon;
  \param[in]      &udiag    Address of vector containing the diagonal of the hessian at the old iterate.
 */
template<class V, class M, class FUNC> void Uncmin<V, M, FUNC>::tregup(V &x, double &f, V &g, M &a,
    V &sc, V &sx, bool &nwtake, double &dlt, int &iretcd, V &xplsp, double &fplsp, V &xpls,
    double &fpls, bool &mxtake, int method, V &udiag)
{
  int i, j;
  double rln, temp, dltf, slp, dltmp, dltfp;

  mxtake = false;

  for (i = 1; i <= n; i++)
  {
    xpls(i) = x(i) + sc(i);
  }
  fpls = minclass->f_to_minimize(xpls);

  dltf = fpls - f;

  slp = ddot(g, sc); // Blas_f77.ddot_f77(n,g,1,sc,1);

  if (iretcd == 4)
    fplsp = 0.0;

  if ((iretcd == 3) && ((fpls >= fplsp) || (dltf > .0001*slp)))
  {
    // RESET XPLS TO XPLSP AND TERMINATE GLOBAL STEP
    iretcd = 0;

    for (i = 1; i <= n; i++)
    {
      xpls(i) = xplsp(i);
    }
    fpls = fplsp;
    dlt *= .5;
  }
  else
  {
    // FPLS TOO LARGE
    if (dltf > .0001*slp)
    {
      rln = 0.0;

      for (i = 1; i <= n; i++)
      {
        rln = max(rln, std::fabs(sc(i))/ max(std::fabs(xpls(i)), 1.0/sx(i)));
      }

      if (rln < steptl)
      {
        // CANNOT FIND SATISFACTORY XPLS SUFFICIENTLY DISTINCT FROM X
        iretcd = 1;
      }
      else
      {
        // REDUCE TRUST REGION AND CONTINUE GLOBAL STEP
        iretcd = 2;

        dltmp = -slp*dlt/(2.0*(dltf - slp));

        if (dltmp < .1*dlt)
        {
          dlt *= .1;
        }
        else
        {
          dlt = dltmp;
        }
      }
    }
    else
    {
      // FPLS SUFFICIENTLY SMALL
      dltfp = 0.0;

      if (method == 2)
      {
        for (i = 1; i <= n; i++)
        {
          temp = 0.0;

          for (j = i; j <= n; j++)
          {
            temp += a(j, i)*sc(j);
          }
          dltfp += temp*temp;
        }
      }
      else
      {
        for (i = 1; i <= n; i++)
        {
          dltfp += udiag(i)*sc(i)*sc(i);
          temp = 0.0;

          for (j = i+1; j <= n; j++)
          {
            temp += a(i, j)*sc(i)*sc(j);
          }
          dltfp += 2.0*temp;
        }
      }
      dltfp = slp + dltfp/2.0;

      if ((iretcd != 2) && (std::fabs(dltfp - dltf) <= .1*std::fabs(dltf)) && (!nwtake) && (dlt
          <= .99*stepmx))
      {
        // DOUBLE TRUST REGION AND CONTINUE GLOBAL STEP
        iretcd = 3;

        for (i = 1; i <= n; i++)
        {
          xplsp(i) = xpls(i);
        }
        fplsp = fpls;
        dlt = min(2.0*dlt, stepmx);
      }
      else
      {
        // ACCEPT XPLS AS NEXT ITERATE.  CHOOSE NEW TRUST REGION.
        iretcd = 0;

        if (dlt > .99*stepmx)
          mxtake = true;

        if (dltf >= .1*dltfp)
        {
          // DECREASE TRUST REGION FOR NEXT ITERATION
          dlt *= .5;
        }
        else
        {
          // CHECK WHETHER TO INCREASE TRUST REGION FOR NEXT ITERATION
          if (dltf <= .75*dltfp)
            dlt = min(2.0*dlt, stepmx);
        }
      }
    }
  }
}

/*! 
  \brief
  Calculate the dot product of two vectors. 
  
  \section template_args Template Parameters
  
  \param M  SCPPNT Matrix type.
  \param V  SCPPNT Vector type.
  \param FUNC Type of function to minimize.

  \section function_args Function Parameters
   
  \param[in]  &x  Address of first vector to multiply.
  \param[in]  &y  Address of second vector to multiply.
  */
template<class V, class M, class FUNC> double Uncmin<V, M, FUNC>::ddot(V &x, V &y)
{
  typename V::iterator ix = x.begin();
  typename V::iterator iy = y.begin();
  double prod = 0.0;

  // assumes x and y are the same size
  for (int n = x.size(); n--; ++ix, ++iy)
  {
    prod += *ix * *iy;
  }

  return prod;
}

/*! 
  \brief
  Calculate the Euclidean norm of a vector. 
  
  It is a translation from FORTRAN to Java of the LINPACK function
  DNRM2.

  In the LINPACK listing DNRM2 is attributed to C.L. Lawson
  with a date of January 8, 1978.  The routine below is based
  on a more recent DNRM2 version that is attributed in LAPACK
  documentation to Sven Hammarling.

  Translated by Steve Verrill, February 25, 1997.

  \section template_args Template Parameters
  
  \param M  SCPPNT Matrix type.
  \param V  SCPPNT Vector type.
  \param FUNC Type of function to minimize.
  
  \section function_args Function Parameters
   
  \param[in]  &x  Address of first vector to multiply.
 */
template<class V, class M, class FUNC> double Uncmin<V, M, FUNC>::dnrm2(V &x)
{

  double absxi, norm, scale, ssq, fac;
  int ix;
  int nelem = x.size();

  if (nelem < 1)
  {
    norm = 0.0;
  }
  else if (nelem == 1)
  {
    norm = std::fabs(x(1));
  }
  else
  {
    scale = 0.0;
    ssq = 1.0;

    for (ix = 1; ix <= nelem; ++ix)
    {
      if (x(ix) != 0.0)
      {
        absxi = std::fabs(x(ix));

        if (scale < absxi)
        {
          fac = scale/absxi;
          ssq = 1.0 + ssq*fac*fac;
          scale = absxi;
        }
        else
        {
          fac = absxi/scale;
          ssq += fac*fac;
        }
      }
    }
    norm = scale*std::sqrt(ssq);
  }
  return norm;
}

/*! 
  \brief
  Maximum of two double precision scalars.
  
  \section template_args Template Parameters
  
  \param M  SCPPNT Matrix type.
  \param V  SCPPNT Vector type.
  \param FUNC Type of function to minimize.
  
  \section function_args Function Parameters
   
  \param[in]  a First scalar to compare.
  \param[in]  b Second scalar to compare.
 */
template<class V, class M, class FUNC> double Uncmin<V, M, FUNC>::max(double a, double b)
{
  return (a > b) ? a : b;
}

/*! 
  \brief
  Minimum of two double precision scalars.
  
  \section template_args Template Parameters
  
  \param M  SCPPNT Matrix type.
  \param V  SCPPNT Vector type.
  \param FUNC Type of function to minimize.
  
  \section function_args Function Parameters
   
  \param[in]  a First scalar to compare.
  \param[in]  b Second scalar to compare.
 */
template<class V, class M, class FUNC> double Uncmin<V, M, FUNC>::min(double a, double b)
{
  return (a > b) ? b : a;
}

#endif // UNCMIN_H_
