/*!\file   scalar_search.hpp
 * \brief  Functions used to perform minimization in 1-dimension
 * \author Renaud Bruneliere
 * \date   15.10.2013
 *
 * The header contains different functions which can be used to find a minimum
 * along a line like the Golden section search method, a method very similar to
 * the bisection method used in root finding.
 */

#ifndef SCALAR_SEARCH_HPP
#define SCALAR_SEARCH_HPP

#include <utility>
#include "optimize/function.hpp"

namespace optimize {

//! Bracket the minimum point along a line.
/*!
  \param func C-like function pointer or boost::function (to minimize)
  \param x0   starting value
  \param h    starting step size
  \param a    lower bound
  \param b    upper bound
  \param c    step size expansion coefficient
  \param max_iter maximum number of iterations
  \return 0 if a minimum has been found
 */
// @{
int bracket(const t_func & func, double x0, double h,
            double & a, double & b,
            double c=1.6118033989, unsigned int max_iter= 100);
int bracket(t_boofunc & func, double x0, double h,
            double & a, double & b,
            double c=1.6118033989, unsigned int max_iter= 100);
// @}

//! Bracket the minimum point along a line.
/*!
  \param func static wrapper function to callback a class method
  \param obj  object pointer
  \param x0   starting value
  \param h    starting step size
  \param a    lower bound
  \param b    upper bound
  \param c    step size expansion coefficient
  \param max_iter maximum number of iterations
  \return 0 if a minimum has been found
 */
int bracket(const t_func_obj & func, void * obj, 
            double x0, double h,
            double & a, double & b,
            double c=1.6118033989, unsigned int max_iter= 100);

//! Bracket the minimum point along a line.
/*!
  \param func pointer to the objective function
  \param x0   starting value
  \param h    starting step size
  \param a    lower bound
  \param b    upper bound
  \param c    step size expansion coefficient
  \param max_iter maximum number of iterations
  \return 0 if a minimum has been found
 */
int bracket(function * const func, double x0, double h, 
            double & a, double & b,
            double c=1.6118033989, unsigned int max_iter= 100);

//! Find a minimum along a line with the golden section search method.
/*!
  \param func C-like function pointer or boost::function (to minimize)
  \param a    lower bound
  \param b    upper bound
  \param fmin value of the function at minimum
  \param tol  tolerance error
  \return coordinate of the local minimum
 */
// @{
double minimize_1d_golden(const t_func & func, double a, double b,
                          double & fmin, double tol=1.e-9);
double minimize_1d_golden(t_boofunc & func, double a, double b,
                          double & fmin, double tol=1.e-9);
// @}

//! Find a minimum along a line with the golden section search method.
/*!
  \param func static wrapper function to callback a class method
  \param obj  object pointer
  \param a    lower bound
  \param b    upper bound
  \param fmin value of the function at minimum
  \param tol  tolerance error
  \return coordinate of the local minimum
 */
double minimize_1d_golden(const t_func_obj & func, void * obj, 
                          double a, double b,
                          double & fmin, double tol=1.e-9);

//! Find a minimum along a line with the golden section search method.
/*!
  \param func pointer to the objective function
  \param a    lower bound
  \param b    upper bound
  \param fmin value of the function at minimum
  \param tol  tolerance error
  \return coordinate of the local minimum
 */
double minimize_1d_golden(function * const func,
                          double a, double b,
                          double & fmin, double tol=1.e-9);

} // namespace optimize

#endif // SCALAR_SEARCH_HPP
