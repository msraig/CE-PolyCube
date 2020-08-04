/*! \file Fmin.h
 
  \brief
  Univariate function minimization by Brent's method.
  Release 080309.

	A translation of Fmin.java by Steve Verrill to C++ 
	Version 0.000318 (March 18, 2000)
	by Brad Hanson (http://www.b-a-h.com/index.html)
	
	This software is available at 
	http://www.smallwaters.com/software/cpp/uncmin.html
	
	Fmin.java can be obtained from http://www1.fpl.fs.fed.us/optimization.html
	Information about Fmin.java is given below.
	
    Fmin.java copyright claim:
    This software is based on the public domain fmin routine.
    The FORTRAN version can be found at
    www.netlib.org
    This software was translated from the FORTRAN version
    to Java by a US government employee on official time.  
    Thus this software is also in the public domain.
    The translator's mail address is:
    Steve Verrill 
    USDA Forest Products Laboratory
    1 Gifford Pinchot Drive
    Madison, Wisconsin
    53705
    The translator's e-mail address is:
    steve@ws13.fpl.fs.fed.us
    
    (c) 2008, Werner Wothke, Fmin documentation.
    
  ***********************************************************************
  DISCLAIMER OF WARRANTIES:
  THIS SOFTWARE IS PROVIDED "AS IS" WITHOUT WARRANTY OF ANY KIND. 
  THE TRANSLATOR DOES NOT WARRANT, GUARANTEE OR MAKE ANY REPRESENTATIONS 
  REGARDING THE SOFTWARE OR DOCUMENTATION IN TERMS OF THEIR CORRECTNESS, 
  RELIABILITY, CURRENTNESS, OR OTHERWISE. THE ENTIRE RISK AS TO 
  THE RESULTS AND PERFORMANCE OF THE SOFTWARE IS ASSUMED BY YOU. 
  IN NO CASE WILL ANY PARTY INVOLVED WITH THE CREATION OR DISTRIBUTION 
  OF THE SOFTWARE BE LIABLE FOR ANY DAMAGE THAT MAY RESULT FROM THE USE 
  OF THIS SOFTWARE.
  Sorry about that.
  ***********************************************************************
  History:
  Date        Translator        Changes
  3/09/08     Werner Wothke     Doxygen-compatible documentation.
  
  3/24/98     Steve Verrill     Translated
   */
  /**
   *
   *<p>
   *This class was translated by a statistician from the FORTRAN 
   *version of fmin.  It is NOT an official translation.  When 
   *public domain Java optimization routines become available 
   *from professional numerical analysts, then <b>THE CODE PRODUCED
   *BY THE NUMERICAL ANALYSTS SHOULD BE USED</b>.
   *
   *<p>
   *Meanwhile, if you have suggestions for improving this
   *code, please contact Steve Verrill at steve@ws13.fpl.fs.fed.us.
   *
   *@author Steve Verrill
   *@version .5 --- March 24, 1998
   * 
   */
#include <limits>
#include <cmath>
/*!
  \brief 
  This function performs a uni-dimensional minimization by Brent's method.
  
  Brent's method combines a golden-section search and parabolic interpolation.  
  The following introductory comments are copied from the FORTRAN version.
  
  Convergence is never much slower than that for a Fibonacci search.  If f 
  has a continuous second derivative which is positive at the minimum (which 
  is not at ax or bx), then convergence is superlinear, and usually of the 
  order of about 1.324.
  
  The function f is never evaluated at two points closer together
  than eps*abs(fmin)+(tol/3), where eps is approximately the square
  root of the relative machine precision.  If f is a unimodal
  function and the computed values of f are always unimodal when
  separated by at least eps*abs(x)+(tol/3), then fmin approximates
  the abcissa of the global minimum of f on the interval (ax,bx) with
  an error less than 3*eps*abs(fmin)+tol.  If f is not unimodal,
  then fmin may approximate a local, but perhaps non-global, minimum to
  the same accuracy.
  
  This function subprogram is a slightly modified version of the
  Algol 60 procedure localmin given in Richard Brent, Algorithms For
  Minimization Without Derivatives, Prentice-Hall, Inc. (1973).
  
   This method is a translation from FORTRAN to Java of the Netlib 
   function fmin.  In the Netlib listing no author is given.
   
   Translated by Steve Verrill, March 24, 1998.
   
   Translated from Java to C++ by Brad Hanson (bradh@pobox.com)
        
  \section template_args Template Arguments
   
  \param T    Floating point type (e.g., double, float)
  \param FUNC Function object type giving function to minimize.

  \section func_args Function Arguments
   
  \param[in]	a 	Left endpoint of initial interval
  \param[in]	b 	Right endpoint of initial interval
  \param[in]	f		A class that overloads the function call operator
 			[operator()()] to implement the function f(x) for any x
      in the interval (a,b)
  \param[in]	tol	Desired length of the interval in which
 			the minimum will be determined to lie
   		(This should be greater than, roughly, 3.0e-8.)
 */
template<class T, class FUNC>
T Fmin(T a, T b, FUNC &f, T tol) 
{
      T c,d,e,eps,xm,p,q,r,tol1,t2,
             u,v,w,fu,fv,fw,fx,x,tol3;
      c = .5*(3.0 - std::sqrt(5.0));
      d = 0.0;
      // 1.1102e-16 is machine precision
      eps = std::numeric_limits<T>::epsilon();
      tol1 = eps + 1.0;
      eps = std::sqrt(eps);
      v = a + c*(b-a);
      w = v;
      x = v;
      e = 0.0;
      fx = f(x);
      fv = fx;
      fw = fx;
      tol3 = tol/3.0;
      xm = .5*(a + b);
      tol1 = eps*std::fabs(x) + tol3;
      t2 = 2.0*tol1;
      // main loop
      while (std::fabs(x-xm) > (t2 - .5*(b-a))) {
         p = q = r = 0.0;
         if (std::fabs(e) > tol1) {
           // fit the parabola
            r = (x-w)*(fx-fv);
            q = (x-v)*(fx-fw);
            p = (x-v)*q - (x-w)*r;
            q = 2.0*(q-r);
            if (q > 0.0) {
               p = -p;
            } else {
               q = -q;
            }
            r = e;
            e = d;         
            // brace below corresponds to statement 50
         }
         if ((std::fabs(p) < std::fabs(.5*q*r)) &&
             (p > q*(a-x)) &&
             (p < q*(b-x))) {
           // a parabolic interpolation step
            d = p/q;
            u = x+d;
            // f must not be evaluated too close to a or b
            if (((u-a) < t2) || ((b-u) < t2)) {
               d = tol1;
               if (x >= xm) d = -d;
            }
            // brace below corresponds to statement 60
         } else {
           // a golden-section step
            if (x < xm) {
               e = b-x;
            } else {
               e = a-x;
            }
            d = c*e;
         }
         // f must not be evaluated too close to x
         if (std::fabs(d) >= tol1) {
            u = x+d;
         } else {
            if (d > 0.0) {
               u = x + tol1;
            } else {
               u = x - tol1;
            }
         }
         fu = f(u);
         // Update a, b, v, w, and x
         if (fx <= fu) {
            if (u < x) {
               a = u;
            } else {
               b = u;
            }
            // brace below corresponds to statement 140
         }
         if (fu <= fx) {
            if (u < x) {
               b = x;
            } else {
               a = x;
            }
            v = w;
            fv = fw;
            w = x;
            fw = fx;
            x = u;
            fx = fu;
            xm = .5*(a + b);
            tol1 = eps*std::fabs(x) + tol3;
            t2 = 2.0*tol1;
            // brace below corresponds to statement 170
         } else {
            if ((fu <= fw) || (w == x)) {
               v = w;
               fv = fw;
               w = u;
               fw = fu;
               xm = .5*(a + b);
               tol1 = eps*std::fabs(x) + tol3;
               t2 = 2.0*tol1;
            } else if ((fu > fv) && (v != x) && (v != w)) {
               xm = .5*(a + b);
               tol1 = eps*std::fabs(x) + tol3;
               t2 = 2.0*tol1;
            } else {
               v = u;
               fv = fu;
               xm = .5*(a + b);
               tol1 = eps*std::fabs(x) + tol3;
               t2 = 2.0*tol1;
            }
         }
         // brace below corresponds to statement 190
      }
      return x;
}   
                                    
